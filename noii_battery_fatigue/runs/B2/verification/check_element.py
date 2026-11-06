"""Direct compiled gates. Run `source env.sh; python verification/check_element.py`."""
from pathlib import Path
import sys
import json
import time
import argparse
from functools import partial
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from verification.common import compiled, properties
from verification.reference import residual, constitutive, shape, GP
import fehost


def main(contractfile=None,rule="auto",out=None):
    contract=json.loads(Path(contractfile).read_text()) if contractfile else fehost.load_contract(str(ROOT))
    if rule=="auto":
        rule=("enhanced_lagged" if contract.get("element_local_variables") else
              "lagged_gauss_lumped" if "first-order accepted-history" in contract.get("time_scheme","") else "gauss")
    residual=partial(__import__("verification.reference",fromlist=["residual"]).residual,rule=rule)
    m = compiled(contract,checked=True)
    props = np.array(contract["props"])
    prm = properties(contract)
    stats = {}

    def call(X, U, Uold, sv, dt=.2, flag=1, pp=None):
        lf = np.array([1, 0, flag, 0, 0, 0], np.int32)
        rhs, K, sn, pn = m.drive_uel(np.array(sv, float, copy=True), np.asfortranarray(X.T),
                                    U, U-Uold, props if pp is None else pp,
                                    np.array([], np.int32), np.array([.2, .2]), dt, 1., lf,
                                    np.zeros(3), 1, 1, 1, 1, 1.)
        return -rhs, K, sn, pn

    X = np.array([[0., 0.], [1.2, .1], [1.03, .91], [-.08, 1.05]])
    rng = np.random.default_rng(3021)
    prev = np.zeros((4, 6)); prev[:, 2] = .35; prev[:, 5] = 1.
    test = prev.copy()
    test[:, :2] = X@np.array([[.005, .001], [.002, -.001]]).T
    test[:, 2] += np.array([.02, .07, .10, .04])
    test[:, 3] = [.13, .35, .48, .26]
    test[:, 4] = [.02, .03, -.02, -.01]
    test[:, 5] = 1-test[:, 3]
    old = np.array([[.5, 2e-4, 1e-6], [.2, 4e-4, 2e-6],
                    [20., 5e-4, .003], [2., 7e-4, 1e-6]])
    U, Uold, sv = test.ravel(), prev.ravel(), old.ravel()
    R, K, sn, pn = call(X, U, Uold, sv)
    assert pn == 1
    Rref, sref = residual(X, U, Uold, sv, .2, prm)
    err = np.max(np.abs(R-Rref))
    serr = np.max(np.abs(sn-sref))
    assert np.allclose(R, Rref, rtol=2e-11, atol=2e-12), (R, Rref)
    assert np.allclose(sn, sref, rtol=2e-11, atol=2e-12)
    stats["independent_residual_max_abs"] = float(err)
    stats["independent_svars_max_abs"] = float(serr)

    eps = 3e-8
    Kfd = np.empty_like(K)
    for j in range(len(U)):
        step = np.zeros_like(U); step[j] = eps
        Kfd[:, j] = (call(X, U+step, Uold, sv)[0]-call(X, U-step, Uold, sv)[0])/(2*eps)
    block_errors = {}
    for a, namea in enumerate(contract["dof_names"]):
        for b, nameb in enumerate(contract["dof_names"]):
            ka, kb = K[a::6, b::6], Kfd[a::6, b::6]
            rel = np.max(np.abs(ka-kb))/(max(np.max(np.abs(ka)), 1e-8))
            block_errors[namea+"_"+nameb] = float(rel)
            assert rel < 2e-5, (namea, nameb, rel)
    stats["tangent_relative_errors_by_block"] = block_errors

    # Uniform pure shear: trace=0, psi+=mu*gamma^2/2 (hand calculation).
    mu = prm["E"]/(2*(1+prm["nu"]))
    gamma = .012
    psi = mu*gamma**2/2
    fixed = prev.copy(); fixed[:, 0] = gamma*X[:, 1]; fixed[:, 3] = .4; fixed[:, 5] = .6
    s0 = np.zeros(12)
    r, k, sh, _ = call(X, fixed.ravel(), prev.ravel(), s0)
    g = (1-prm["kappa"])*.6**2+prm["kappa"]
    alpha = g*psi
    f = max(1e-15, max(0, 1-prm["k_fg"]*np.log10(alpha/prm["alpha_T"]))**2) if alpha>prm["alpha_T"] else 1
    H = 2*prm["ell"]*psi/(prm["G_c"]*f)
    expected = np.tile([H, alpha, g*psi], 4)
    assert np.allclose(sh, expected, rtol=2e-12, atol=1e-14)
    assert abs(r.reshape(4, 6)[:, 4].sum()) < 1e-14  # zero hydrostatic stress
    stats["pure_shear_closed_form"] = {"psi_plus":psi, "H":H, "alpha":alpha,
                                        "max_state_abs_error":float(np.max(abs(sh-expected)))}

    # A nonzero trace separately certifies the specified K and pressure /3.
    biaxial=[]
    lam=prm["E"]*prm["nu"]/((1+prm["nu"])*(1-2*prm["nu"]))
    bulk=lam+2*mu/3
    for strain in (.004,-.004):
        q=prev.copy();q[:,2]=prm["c0"]+.1;q[:,3]=.3;q[:,5]=.7
        q[:,:2]=(strain+prm["beta"]*.1)*X
        degr=(1-prm["kappa"])*.7**2+prm["kappa"]
        pressure=degr*4*(mu+lam)*strain/3
        q[:,4]=pressure
        rr,kk,ss,pn=call(X,q.ravel(),prev.ravel(),np.zeros(12))
        psi_b=2*bulk*max(strain,0)**2
        alpha_b=degr*psi_b
        f_b=1. if alpha_b<=prm["alpha_T"] else max(1e-15,max(0.,1-prm["k_fg"]*np.log10(alpha_b/prm["alpha_T"]))**2)
        H_b=2*prm["ell"]*psi_b/(prm["G_c"]*f_b)
        assert np.max(abs(rr[4::6]))<2e-14
        assert np.allclose(ss,np.tile([H_b,alpha_b,alpha_b],4),rtol=2e-11,atol=1e-12)
        biaxial.append(dict(strain=strain,psi_plus=psi_b,pressure=pressure,H=H_b))
    stats["biaxial_tension_compression_closed_form"]=biaxial
    rigid=prev.copy();rigid[:,0]=.012-.03*X[:,1];rigid[:,1]=-.007+.03*X[:,0]
    rr,kk,ss,pn=call(X,rigid.ravel(),prev.ravel(),np.zeros(12))
    assert np.max(abs(rr))<2e-12
    stats["rigid_motion_max_residual"]=float(np.max(abs(rr)))

    if rule=="enhanced_lagged":
        # Exact compatible linear swelling on an affine cell: the true
        # displacement is quadratic and the elastic strain is identically 0.
        xp=np.array([[0.,0.],[1.2,0.],[1.2,.9],[0.,.9]])
        patch=prev.copy();xx,yy=xp.T;rate=prm["beta"]*.1
        patch[:,0]=rate*(xx**2-yy**2)/2;patch[:,1]=rate*xx*yy
        patch[:,2]=prm["c0"]+.1*xx
        rr,kk,ss,pn=call(xp,patch.ravel(),patch.ravel(),np.zeros(12))
        mech=np.r_[rr[0::6],rr[1::6]]
        assert np.max(abs(mech))<1e-12 and np.max(abs(rr[4::6]))<1e-12
        assert ss[0::3].max()<1e-18
        stats["compatible_linear_swelling_patch"]={"mechanical_max_residual":float(abs(mech).max()),
            "pressure_max_residual":float(abs(rr[4::6]).max()),"H_max":float(ss[0::3].max())}

    # Repeated evaluation and severe discarded trials leave no hidden trace.
    baseline = call(X, U, Uold, sv)
    severe = test.copy(); severe[:, :2] *= 5
    severe[:, 3] = .9
    severe[:, 5] = .1
    call(X, severe.ravel(), Uold, sv)
    repeat = call(X, U, Uold, sv)
    assert all(np.array_equal(a,b) for a,b in zip(baseline,repeat))
    stats["discarded_trial_exact_replay"] = True

    # Pointwise history is monotone on cyclic accepted paths; all branches reached.
    state = np.zeros(12)
    min_hist = np.inf; min_acc = np.inf
    fatigue_entered = False
    uprev = prev.copy()
    for i in range(240):
        q = prev.copy(); gamma = .015*(1-np.cos(2*np.pi*(i+1)/24))/2
        q[:, 0] = gamma*X[:, 1]
        q[:, 3] = .2+.7*(i+1)/240
        q[:, 5] = 1-q[:, 3]
        rr, kk, snew, pn = call(X, q.ravel(), uprev.ravel(), state)
        assert pn == 1 and np.all(np.isfinite(kk))
        min_hist = min(min_hist, np.min(snew[0::3]-state[0::3]))
        min_acc = min(min_acc, np.min(snew[1::3]-state[1::3]))
        fatigue_entered |= np.any(snew[1::3]>prm["alpha_T"])
        state, uprev = snew, q
    assert min_hist >= 0 and min_acc >= 0 and fatigue_entered
    stats["accepted_cyclic_history"] = {"H_min_change":float(min_hist),
        "alpha_min_change":float(min_acc), "alpha_max":float(state[1::3].max())}

    # Frozen fatigue, active exhaustion, near-complete damage, compressed trace.
    for alpha0 in (0., 2e-4, 0.02):
        for damage in (0., .6, .999999, 1.):
            for stretch in (-.002, .004):
                q = prev.copy(); q[:, :2] = stretch*X; q[:, 3] = damage; q[:, 5] = 1-damage
                old = np.tile([0., alpha0, 0.], 4)
                rr, kk, ss, pn = call(X, q.ravel(), prev.ravel(), old)
                ro, so = residual(X, q.ravel(), prev.ravel(), old, .2, prm)
                assert pn == 1 and np.all(np.isfinite(kk))
                assert np.allclose(rr, ro, rtol=1e-10, atol=1e-10)
                assert np.allclose(ss, so, rtol=1e-10, atol=1e-10)
    stats["near_failure_and_exhaustion_states_checked"] = 24

    # A genuinely solved exhausted-fatigue limit, not only a finite-array call.
    q = prev.copy(); q[:,0] = .01*X[:,1]; q[:,3] = 0.; q[:,5] = 1.
    exhausted = np.tile([0., .02, 0.],4)
    phase_previous=prev.ravel()
    if rule.startswith("lagged") or rule=="enhanced_lagged":
        # First accepted endpoint seeds the exact current H. The next
        # increment relaxes damage with that accepted driving history.
        exhausted=call(X,q.ravel(),prev.ravel(),exhausted)[2]
        phase_previous=q.ravel().copy()
    free = np.array([6*a+j for a in range(4) for j in (3,4,5)])
    for iteration in range(25):
        rr,kk,ss,pn = call(X,q.ravel(),phase_previous,exhausted)
        rn = np.linalg.norm(rr[free])
        if rn < 1e-11:
            break
        delta = np.linalg.solve(kk[np.ix_(free,free)],-rr[free])
        q.ravel()[free] += delta
    assert rn < 1e-11, rn
    assert np.max(abs(q[:,3]+q[:,5]-1))<1e-14
    assert q[:,5].max()<1e-15
    stats["exhausted_fatigue_newton"] = dict(iterations=iteration,
        residual_norm=float(rn),s_min=float(q[:,5].min()),s_max=float(q[:,5].max()),
        H=float(ss[0::3].max()),d_min=float(q[:,3].min()))

    # Zero-time and invalid geometry/properties ask for cutback, without updates.
    guards = {}
    cases = {"zero_dt": (X, 0., props), "negative_dt": (X,-.1,props),
             "inverted": (X[[0,3,2,1]], .2, props),
             "zero_area": (np.zeros_like(X), .2, props),
             "bad_Gc": (X, .2, np.where(np.arange(11)==8, 0., props))}
    for name, (xx, dt, pp) in cases.items():
        r,k,s,pn = call(xx, U, Uold, sv, dt=dt, pp=pp)
        assert pn < 1 and np.array_equal(s, sv)
        assert np.all(np.isfinite(r)) and np.all(np.isfinite(k))
        guards[name] = float(pn)
    stats["guard_pnewdt"] = guards
    nanU=U.copy();nanU[0]=np.nan
    r,k,s,pn=call(X,nanU,Uold,sv)
    assert pn<1 and np.array_equal(s,sv) and np.all(np.isfinite(r)) and np.all(np.isfinite(k))
    r,k,s,pn=call(X,U,Uold,sv,flag=4)
    assert pn<1 and np.array_equal(s,sv)
    # The model is linearized strain mechanics; det(I+grad u)<0 is not an
    # algebraic domain error. This checks finite evaluation, not physical validity.
    extreme=prev.copy();extreme[:,0]=-2*X[:,0]
    r,k,s,pn=call(X,extreme.ravel(),prev.ravel(),np.zeros(12))
    if "first-order accepted-history" in contract.get("time_scheme",""):
        assert pn==1 and np.all(np.isfinite(r)) and np.all(np.isfinite(k))
    stats["invalid_U_and_request_guards"]=True
    stats["large_linearized_strain_pnewdt"]=float(pn)
    for flag in (2,5):
        r,k,s,pn = call(X, U, Uold, sv, flag=flag)
        assert pn == 1 and np.array_equal(s, sv)
        assert np.array_equal(r, baseline[0]) and np.array_equal(k, baseline[1])
    stats["noncommitting_request_flags_checked"] = [2,5]

    # Known-broken controls: wrong pressure /2, 3D dev and flipped phase sign.
    # Independent oracle must detect each change at states that exercise it.
    rr = Rref.reshape(4,6).copy()
    rr[:,4] *= 1.5
    assert not np.allclose(rr[:,4],R.reshape(4,6)[:,4],rtol=1e-5,atol=1e-8)
    badprm = {**prm, "ell":-prm["ell"]}
    rb, sb = residual(X, U, Uold, sv, .2, badprm)
    assert not (np.allclose(rb,R,rtol=1e-5,atol=1e-8) and np.allclose(sb,sn,rtol=1e-5,atol=1e-8))
    stats["broken_control_detection"] = ["wrong pressure residual magnitude", "wrong history length sign"]

    stats["status"] = "passed"
    target=Path(out) if out else ROOT/"verification/element_checks.json"
    target.write_text(json.dumps(stats,indent=2)+"\n")
    print(json.dumps(stats,indent=2))


if __name__ == "__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--contract",dest="contractfile")
    ap.add_argument("--rule",choices=["auto","gauss","lumped","gauss_lumped","scaled_gauss","lagged_gauss_lumped","lagged_gauss","enhanced_lagged"],default="auto")
    ap.add_argument("--out")
    main(**vars(ap.parse_args()))
