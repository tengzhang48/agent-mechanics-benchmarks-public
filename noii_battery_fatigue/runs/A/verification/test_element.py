"""Finite-difference, analytic patch, conservation and statelessness tests.

Run with the pinned environment, TMPDIR=$PWD/tmp and bytecode disabled.
"""
import json
import numpy as np
from common import ROOT, load_element, raw_call, write_json


def tangent_check(el, x, u, old, sv):
    r, k, _, pn = raw_call(el, x, u, old, sv)
    assert pn >= 1
    errors = []
    for eps in [1e-5, 1e-6, 1e-7]:
        fd = np.empty_like(k)
        for j in range(len(u)):
            inc = np.zeros_like(u)
            inc[j] = eps
            rp, _, _, pp = raw_call(el, x, u+inc, old, sv)
            rm, _, _, pm = raw_call(el, x, u-inc, old, sv)
            assert min(pp, pm) >= 1
            fd[:, j] = (rp-rm)/(2*eps)
        errors.append(float(np.linalg.norm(fd-k)/np.linalg.norm(k)))
    assert min(errors) < 2e-8, errors
    return errors


def main():
    el, contract = load_element()
    x = np.array([[0., 0.], [1.2, .08], [1.08, .93], [-.13, 1.04]])
    rng = np.random.default_rng(821)
    uv = np.zeros((4, 6))
    uv[:, :2] = x @ np.array([[.025, .008], [-.009, -.006]])
    uv[:, :2] += rng.normal(0, .001, (4, 2))
    uv[:, 2] = [.43, .51, .54, .46]
    uv[:, 3] = [.15, .22, .32, .23]
    uv[:, 4] = [.007, -.01, .02, .003]
    uv[:, 5] = [.17, .3, .4, .22]
    u = uv.ravel()
    old = u.copy()
    old[2::6] -= .03
    active_errors = tangent_check(el, x, u, old, np.zeros(27))
    sv = np.tile([1e6, 1e-6, .1], 9)
    frozen_errors = tangent_check(el, x, u, old, sv)
    _,_,frozen_sn,_=raw_call(el,x,u,old,sv)
    assert np.array_equal(frozen_sn[0::3],sv[0::3])
    assert np.array_equal(frozen_sn[1::3],sv[1::3])
    floor_errors=tangent_check(el,x,u,old,np.tile([0.,.2,.1],9))

    # Uniform physical strain is exactly reproduced even on nonaffine quads.
    p = dict(zip(contract['prop_names'], contract['props']))
    mu = p['E']/(2*(1+p['nu']))
    lam = p['E']*p['nu']/((1+p['nu'])*(1-2*p['nu']))
    bulk = lam + 2*mu/3
    strain = np.array([[.003, .0015], [.0015, -.0004]])
    sw = p['beta']*(.46-p['c0'])
    ee = strain-sw*np.eye(2)
    tr = np.trace(ee)
    stress0 = 2*mu*ee+lam*tr*np.eye(2)
    dev = ee-.5*tr*np.eye(2)
    psi = .5*bulk*max(tr, 0)**2+mu*np.sum(dev*dev)
    h = 2*p['ell']*psi/p['G_c']
    d = (1-p['kappa'])*h/(1+(1-p['kappa'])*h)
    g = (1-p['kappa'])*(1-d)**2+p['kappa']
    patch = np.column_stack([x @ strain.T + [.02, -.03],
                             np.full(4, .46), np.full(4, d),
                             np.full(4, g*np.trace(stress0)/3), np.full(4, d)]).ravel()
    # Large alpha_T switches off fatigue, so the constant patch is analytic.
    el.props[9] = 1e3
    patch_sv=np.zeros(27); patch_sv[0::3]=h
    r, _, sn, pn = raw_call(el, x, patch, patch,patch_sv)
    assert pn >= 1
    assert np.max(np.abs(r.reshape(4,6)[:,2:])) < 1e-12
    np.testing.assert_allclose(sn[0::3], h, rtol=2e-13, atol=1e-13)
    np.testing.assert_allclose(sn[1::3], g*psi, rtol=2e-13, atol=1e-15)
    np.testing.assert_allclose(sn[2::3], g*psi, rtol=2e-13, atol=1e-15)
    # Independent boundary integration of constant stress * outward normal.
    forces = np.zeros((4,2))
    for a in range(4):
        b = (a+1)%4
        edge = x[b]-x[a]
        normal_ds = np.array([edge[1], -edge[0]])
        forces[a] += .5*g*stress0 @ normal_ds
        forces[b] += .5*g*stress0 @ normal_ds
    patch_error = float(np.max(np.abs(r.reshape(4,6)[:,:2]-forces)))
    assert patch_error < 1e-13

    # Equibiaxial tension isolates the specified K and positive-trace split.
    tensile=patch.reshape(4,6).copy()
    tensile[:,:2]=.01*x
    tensile[:,2]=p['c0']
    tensile_psi=.5*bulk*.02**2
    tensile_h=2*p['ell']*tensile_psi/p['G_c']
    tensile_d=(1-p['kappa'])*tensile_h/(1+(1-p['kappa'])*tensile_h)
    tensile_g=(1-p['kappa'])*(1-tensile_d)**2+p['kappa']
    tensile[:,3]=tensile_d
    tensile[:,4]=tensile_g*2*(mu+lam)*.02/3
    tensile[:,5]=tensile_d
    tensile_sv=np.zeros(27); tensile_sv[0::3]=tensile_h
    rt,_,st,pt=raw_call(el,x,tensile.ravel(),tensile.ravel(),tensile_sv)
    np.testing.assert_allclose(st[0::3],tensile_h,rtol=2e-13)
    np.testing.assert_allclose(st[2::3],tensile_g*tensile_psi,rtol=2e-13)
    assert np.max(np.abs(rt.reshape(4,6)[:,2:]))<1e-12
    floor_state=tensile.copy();floor_state[:,3]=.3;floor_state[:,5]=0.
    floor_sv=np.tile([0.,.2,.1],9)
    el.props[9]=p['alpha_T']
    _,_,sf,pf=raw_call(el,x,floor_state.ravel(),floor_state.ravel(),floor_sv)
    np.testing.assert_allclose(sf[0::3],2*p['ell']*tensile_psi/(p['G_c']*1e-15),rtol=5e-13)
    np.testing.assert_array_equal(sf[1::3],floor_sv[1::3])

    rigid=np.zeros((4,6))
    rigid[:,:2]=np.column_stack([.02-.015*x[:,1],-.03+.015*x[:,0]])
    rigid[:,2]=p['c0']
    rg,_,sg,pg=raw_call(el,x,rigid.ravel(),rigid.ravel())
    assert np.max(np.abs(rg))<1e-13 and np.max(sg)<1e-20

    # Exact endpoint history in a prescribed cyclic homogeneous strain.
    # d is held fixed here to isolate the W/alpha/f/H definitions.
    el.props[9]=p['alpha_T']
    cyc=np.zeros((4,6)); cyc[:,2]=p['c0']; cyc[:,3]=.3
    cyc_sv=np.zeros(27); cyc_old=cyc.ravel().copy()
    amplitude=.004; ncycles=3
    for phase in np.arange(1,16*ncycles+1)*2*np.pi/16:
        strain_amp=amplitude*(1-np.cos(phase))/2
        cyc[:,:2]=strain_amp*x
        _,_,cyc_sv,pncyc=raw_call(el,x,cyc.ravel(),cyc_old,cyc_sv)
        assert pncyc>=1
        cyc_old=cyc.ravel().copy()
    peak_psi=2*bulk*amplitude**2
    peak_g=(1-p['kappa'])*.7**2+p['kappa']
    expected_alpha=ncycles*peak_g*peak_psi
    fatigue=max(1e-15,max(0,1-p['k_fg']*np.log10(expected_alpha/p['alpha_T']))**2)
    expected_h=2*p['ell']*peak_psi/(p['G_c']*fatigue)
    np.testing.assert_allclose(cyc_sv[1::3],expected_alpha,rtol=5e-13)
    np.testing.assert_allclose(cyc_sv[0::3],expected_h,rtol=5e-13)
    assert np.max(np.abs(cyc_sv[2::3]))<1e-25

    el.props[9] = p['alpha_T']
    r, k, sn, pn = raw_call(el, x, u, old)
    # All local flux contributions must sum to zero.
    area = .5*abs(np.dot(x[:,0], np.roll(x[:,1],-1))-
                   np.dot(x[:,1], np.roll(x[:,0],-1)))
    mass_error = abs(r[2::6].sum()-area*.03/el.dt)
    assert mass_error < 1e-13

    # Neither call order nor garbage trial evaluations may alter the answer.
    for _ in range(15):
        trial = u+rng.normal(0,.03,len(u))
        raw_call(el,x,trial,old)
    repeat = raw_call(el, x, u, old)
    assert all(np.array_equal(a,b) for a,b in zip((r,k,sn),repeat[:3]))
    assert np.all(el.svars == 0)
    r0, k0, s0, pn0 = raw_call(el,x,u,old,sv,dt=0)
    assert pn0 < 1 and np.isfinite(r0).all() and np.isfinite(k0).all()
    assert np.array_equal(s0,sv)
    badx=x.copy(); badx[1:]=badx[:0:-1]
    ri, ki, si, pni = raw_call(el,badx,u,old,sv)
    assert pni < 1 and np.array_equal(si,sv)
    assert np.isfinite(ri).all() and np.isfinite(ki).all()
    nonfinite=u.copy();nonfinite[0]=np.nan
    rn,kn,sn_bad,pnn=raw_call(el,x,nonfinite,old,sv)
    assert pnn<1 and np.isfinite(rn).all() and np.isfinite(kn).all() and np.array_equal(sn_bad,sv)

    # The clipped fatigue law permits H > 1/machine epsilon. Source residual
    # stays finite and can represent the exactly saturated state d=z=1.
    sat = uv.copy()
    sat[:,3]=1.; sat[:,5]=1.
    sat[:,4]=0.
    sat = sat.ravel()
    saturated_sv=np.tile([1e22, 1e3, 1e-4],9)
    rs, ks, sns, pns = raw_call(el,x,sat,sat,saturated_sv)
    assert pns >= 1 and np.isfinite(ks).all()
    assert np.max(np.abs(rs.reshape(4,6)[:,[3,5]])) < 1e-16
    assert (sns[0::3] >= saturated_sv[0::3]).all()
    assert (sns[1::3] >= saturated_sv[1::3]).all()
    result = {'active_fatigue_tangent_relative_errors': active_errors,
              'frozen_history_tangent_relative_errors': frozen_errors,
              'fatigue_floor_tangent_relative_errors': floor_errors,
              'affine_patch_force_max_error': patch_error,
              'local_mass_balance_error': float(mass_error),
              'equibiaxial_tensile_energy_and_rigid_motion': True,
              'three_cycle_history_matches_analytic_alpha_H_W': True,
              'fatigue_floor_history_matches_analytic_H': True,
              'repeat_after_trial_calls_bitwise_equal': True,
              'zero_dtime_finite_rejection_unchanged_state': True,
              'reversed_geometry_finite_rejection_unchanged_state': True,
              'nonfinite_input_finite_rejection_unchanged_state': True,
              'H_1e22_damage_saturation_finite': True}
    write_json(ROOT/'verification'/'element_tests.json', result)
    print(json.dumps(result, indent=2))


if __name__=='__main__':
    main()
