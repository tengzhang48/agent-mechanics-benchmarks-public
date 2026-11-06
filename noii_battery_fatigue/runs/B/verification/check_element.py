"""Independent assembled equations, FD Jacobian, invariants, and failure tests.

Usage: python verification/check_element.py
The reference below does not import the generator, its model, or its assembly.
"""
import json
from pathlib import Path

import numpy as np

from common import ROOT, compiled, direct, source_identity, identical_bytes


def quadrature(X):
    corners = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])
    # SVARS integration-point identity is cyclic: --, +-, ++, -+.
    for xi, eta in corners / np.sqrt(3):
        s, t = corners.T
        N = (1+s*xi)*(1+t*eta)/4
        DN = np.array([s*(1+t*eta), t*(1+s*xi)]).T/4
        Jac = X.T @ DN
        B = DN @ np.linalg.inv(Jac)
        yield N, B, np.linalg.det(Jac)


def independent(X, U, previous, sv0, dt, props, *, wrong_pressure=False, wrong_dev=False):
    E, nu, beta, c0, kap, D, chi, ell, Gc, aT, kfg = props[:11]
    theta = props[11]
    mu = E/(2*(1+nu))
    lam = E*nu/((1+nu)*(1-2*nu))
    bulk = lam+2*mu/3
    nodal = U.reshape(4, 6)
    old = previous.reshape(4, 6)
    old_state = sv0.reshape(4, 4)
    R = np.zeros((4, 6))
    sn = np.zeros((4, 4))
    for gp, (N, B, w) in enumerate(quadrature(X)):
        # Matrix definition, intentionally independent of component codegen.
        strain = (nodal[:, :2].T @ B + B.T @ nodal[:, :2])/2
        c, d, q, vfrac = N @ nodal[:, 2:]
        ee = strain-beta*(c-c0)*np.eye(2)
        trace = np.trace(ee)
        effective = 2*mu*ee+lam*trace*np.eye(2)
        g = kap+(1-kap)*vfrac**2
        sigma = g*effective
        pressure = np.trace(effective)/(2 if wrong_pressure else 3)
        if wrong_dev:
            # Deliberately broken 3D deviator with e33=0.
            ee3 = np.zeros((3, 3)); ee3[:2, :2] = ee
            dev = ee3-np.trace(ee3)/3*np.eye(3)
        else:
            dev = ee-trace/2*np.eye(2)
        psi = bulk/2*max(trace, 0)**2+mu*np.sum(dev*dev)
        W = g*psi
        h0, a0, w0, psi_old = old_state[gp]
        alpha = a0+max(0., W-w0)
        f = 1. if alpha <= aT else max(1e-15, max(0., 1-kfg*np.log10(alpha/aT))**2)
        H = max(h0, 2*ell*psi/(Gc*f))
        psi_phase = min(psi, max(psi_old, theta*psi))
        H_phase = max(h0, 2*ell*psi_phase/(Gc*f))
        J = g*D*(nodal[:, 2] @ B)-chi*c*(nodal[:, 4] @ B)
        R[:, :2] += B @ sigma.T*w
        R[:, 2] += (N*(c-N@old[:, 2])/dt+B@J)*w
        R[:, 3] += (ell**2*B@(nodal[:, 3]@B)+N*(nodal[:,3]-(1-kap)*nodal[:,5]*H_phase))*w
        R[:, 4] += N*(q-g*pressure)*w
        R[:, 5] += N*((d-1.)+vfrac)*w
        sn[gp] = H, alpha, W, psi
    return R.ravel(), sn.ravel()


def main():
    module, contract = compiled()
    X = np.array([[.0, .0], [1.1, -.08], [1.04, .91], [-.11, 1.02]])
    rng = np.random.default_rng(4289)
    checks = {}
    props = np.asarray(contract["props"])
    # Fully nonuniform trial, not close to a stationary solution or a symmetry.
    base = np.zeros((4, 6))
    base[:, :2] = X @ np.array([[.018, -.004], [.007, .009]])+rng.normal(0,.001,(4,2))
    base[:, 2] = [.41, .57, .62, .46]
    base[:, 3] = [.21, .34, .4, .25]
    base[:, 4] = [.008, -.004, .019, .03]
    base[:,5] = [.8, .6, .9, .7]
    U = base.ravel()
    prev = U.copy(); prev[2::6] -= [.03, .02, .06, .01]
    sv_cases = {
        "fresh": np.zeros((4,4)),
        "fatigue_loading": np.tile([0., .002, 1e-8, .001], (4,1)),
        "frozen_unloading": np.tile([1e6, .002, .1, .003], (4,1)),
        "fatigue_floor": np.tile([1., .1, 1e-9, .0008], (4,1)),
        "distinct_point_state": np.array([[3., .0002, 1e-7, .0001], [90., .002, .1, .003],
                                           [300., .0003, 2e-6, 0.], [800., .004, .01, .0004]]),
    }
    for label, sv in sv_cases.items():
        rhs, K, out, pn = direct(module, contract, X, U, prev, sv.ravel())
        assert pn == 1.
        rr, ss = independent(X, U, prev, sv.ravel(), .2, props)
        re = np.max(abs(-rhs-rr))/max(1., np.max(abs(rr)))
        se = np.max(abs(out-ss))/max(1., np.max(abs(ss)))
        assert re < 2e-12 and se < 2e-12, (label, re, se)
        Kfd = np.zeros_like(K)
        # Hold the previous solution and state fixed: total DU changes with U.
        h = 2e-7
        for j in range(24):
            v = np.eye(24)[j]*h
            rp = direct(module, contract, X, U+v, prev, sv.ravel())[0]
            rm = direct(module, contract, X, U-v, prev, sv.ravel())[0]
            Kfd[:,j] = -(rp-rm)/(2*h)
        # Row-relative checks keep the diffusion and q rows visible beside huge H.
        row_error = np.max(abs(K-Kfd),axis=1)/np.maximum(np.max(abs(K),axis=1),1e-12)
        assert row_error.max() < 3e-6, (label, row_error)
        checks[label] = dict(residual_relative_error=float(re), state_relative_error=float(se),
                             tangent_max_row_relative_error=float(row_error.max()))

    # Broken controls must visibly fail against independent residual/history.
    correct, correctstate = independent(X,U,prev,np.zeros(16),.2,props)
    wrongp, _ = independent(X,U,prev,np.zeros(16),.2,props,wrong_pressure=True)
    wrongdev, wrongstate = independent(X,U,prev,np.zeros(16),.2,props,wrong_dev=True)
    assert np.max(abs(correct[4::6]-wrongp[4::6])) > 1e-4
    assert np.max(abs(correctstate-wrongstate)) > .01
    checks["broken_controls"] = dict(pressure_divisor_rejected=True, dev3_rejected=True)

    # Homogeneous free swelling, affine patch, and rigid infinitesimal rotation.
    for label, rotation in [("free_swelling",0.),("rigid_rotation",.017)]:
        patch = np.zeros((4,6)); patch[:,2] = .63; patch[:,5] = 1.
        a = props[2]*(.63-props[3])
        patch[:,:2] = a*X + X @ np.array([[0.,rotation],[-rotation,0.]]) + [.13,-.07]
        rhs,K,s,pn = direct(module,contract,X,patch.ravel())
        assert pn == 1 and max(abs(rhs)) < 1e-12 and max(abs(s)) < 1e-20
        checks[label] = dict(residual_inf=float(max(abs(rhs))))

    rhs,K,s,pn = direct(module,contract,X,U,prev)
    mass = sum((N@(U.reshape(4,6)[:,2]-prev.reshape(4,6)[:,2]))*w/.2
               for N,B,w in quadrature(X))
    assert abs(sum(-rhs[2::6])-mass) < 1e-13
    checks["concentration_balance"] = float(abs(sum(-rhs[2::6])-mass))

    # Returning to identical input after arbitrary discarded calls is bitwise exact.
    baseline = direct(module,contract,X,U,prev)
    for _ in range(24):
        trial = U+rng.normal(0,.015,24)
        direct(module,contract,X,trial,prev)
    repeat = direct(module,contract,X,U,prev)
    assert all(identical_bytes(a,b) for a,b in zip(baseline,repeat))
    checks["discarded_trial_hidden_memory"] = "bitwise identical R, K, SVARS and PNEWDT"

    # Zero-time and invalid-state requests are finite, zero, state-preserving rejections.
    oldsv = np.tile([10.,.002,.0001,.0007],4)
    for label, kw in [
        ("zero_time", dict(dt=0.)),
        ("negative_time", dict(dt=-.1)),
        ("inverted_geometry", dict(X=X[[0,3,2,1]])),
        ("zero_Gc", dict(props=np.where(np.arange(len(props))==8,0.,props))),
        ("nonfinite_U", dict(U=np.full(24,np.nan))),
    ]:
        args = dict(X=X,U=U,previous=prev,state=oldsv,dt=.2)
        args.update(kw)
        rr,kk,ss,pp = direct(module,contract,**args)
        assert pp < 1 and np.isfinite(rr).all() and np.isfinite(kk).all()
        assert np.array_equal(ss,oldsv) and not rr.any() and not kk.any()
        checks[label] = dict(pnewdt=float(pp), finite=True, state_preserved=True)
    for request in [2,5]:
        flags=contract["lflags"].copy(); flags[2]=request
        rr,kk,ss,pp=direct(module,contract,X,U,prev,oldsv,flags=flags)
        assert pp==1 and np.array_equal(ss,oldsv)
    checks["noncommitting_lflags"] = [2,5]

    # Exact homogeneous damage solution for a fixed committed H. Report roundoff
    # without disguising it as a proof of global solvability at the fatigue floor.
    limits=[]
    for H in [1.,1e4,1e8,1e12,1e15]:
        ad=(1-props[4])*H
        dstar=ad/(1+ad)
        state=np.tile([H,0.,0.,0.],4)
        x=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        un=np.zeros((4,6));un[:,2]=props[3];un[:,3]=dstar
        un[:,5] = 1./(1.+ad)
        rr,kk,ss,pp=direct(module,contract,x,un.ravel(),state=state)
        limits.append(dict(H=H,dstar=dstar,residual_inf=float(max(abs(rr))),pnewdt=float(pp)))
    checks["large_history_roundoff"] = limits
    checks["element_source_sha256"] = source_identity()
    path=ROOT/"verification/element_checks.json"
    path.write_text(json.dumps(checks,indent=2)+"\n")
    print(json.dumps(checks,indent=2))


if __name__=="__main__":
    main()
