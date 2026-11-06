"""Independent closed-form directional oracle for the lumped psi-lag phase row.

Affine physical strain is exact on the distorted isoparametric element.
Distinct old histories select lagged, predictor-floor, current-energy-cap,
and frozen-H branches at the four points. Expected state, residual, and
directional tangent values use the stated laws only; no model or generated
assembly is imported.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys

import numpy as np

from common import ROOT, direct, fehost


def compiled(root):
    """Compile any workspace element root, independent of its helper scripts."""
    contract = fehost.load_contract(str(root))
    sources = contract['source']
    if isinstance(sources, str):
        sources = [sources]
    digest = hashlib.sha256()
    for name in sources:
        digest.update((root/name).read_bytes())
    digest.update(str(contract.get('compile_flags', [])).encode())
    name = 'verified_uel_'+digest.hexdigest()[:12]
    cached = root/'tmp'/name
    build = cached if list(cached.glob(name+'*.so')) else ROOT/'tmp'/name
    if list(build.glob(name+'*.so')):
        sys.path.insert(0, str(build))
        module = importlib.import_module(name)
    else:
        module = fehost.build_module(str(root), contract, str(build), name)
    return module, contract


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--element-root', type=Path, default=ROOT)
    parser.add_argument('--out', type=Path,
                        default=ROOT/'verification/psilag_phase_oracle.json')
    parser.add_argument('--dt', type=float, default=.2)
    parser.add_argument('--relaxation-form',
                        choices=['accepted_h', 'accepted_reaction', 'current_reaction',
                                 'current_reaction_sqrt'],
                        default='accepted_h',
                        help='Independent mathematical expectation; does not inspect model source')
    args = parser.parse_args()
    assert args.dt > 0.
    root = args.element_root.resolve()
    assert root == ROOT or ROOT in root.parents
    module, contract = compiled(root)
    assert contract['nsvars'] == 16 and contract['dof_names'][-1] == 'vfrac'
    props = np.array(contract['props'], float)
    for name, value in [('G_c', .002), ('alpha_T', .0005)]:
        props[contract['prop_names'].index(name)] = value
    p = dict(zip(contract['prop_names'], props))
    indices = {name: i for i, name in enumerate(contract['dof_names'])}
    dpn = len(indices)
    X = np.array([[0., 0.], [1.1, -.05], [.95, .85], [-.12, 1.03]])
    elastic = np.array([[.03, .012], [.012, -.01]])
    concentration, damage, vfrac = .55, .4, .6
    total = elastic+p['beta']*(concentration-p['c0'])*np.eye(2)
    U = np.zeros((4, dpn))
    U[:, [indices['ux'], indices['uy']]] = X@total.T
    U[:, indices['c']] = concentration
    U[:, indices['d']] = damage
    U[:, indices['vfrac']] = vfrac
    a = 1-p['kappa']
    mu = p['E']/(2*(1+p['nu']))
    lam = p['E']*p['nu']/((1+p['nu'])*(1-2*p['nu']))
    bulk = lam+2*mu/3
    trace = np.trace(elastic)
    psi = .5*bulk*max(trace, 0.)**2+mu*(.5*(elastic[0, 0]-elastic[1, 1])**2+2*elastic[0, 1]**2)
    g = p['kappa']+a*vfrac*vfrac
    W = g*psi
    C = 2*p['ell']/p['G_c']
    theta = p['psi_predictor_floor_fraction']
    relaxation = p.get('phase_relaxation_factor', 0.)
    effective_relaxation = (relaxation/np.sqrt(args.dt)
                            if args.relaxation_form == 'current_reaction_sqrt' else relaxation)
    cap = p.get('phase_relaxation_cap', 1000.)
    U[:, indices['q']] = g*2*(mu+lam)*trace/3
    previous = U.copy()
    previous[:, indices['d']] = [.1, .2, .3, .25]
    previous[:, indices['vfrac']] = 1-previous[:, indices['d']]

    signs = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])
    Q = []
    for xi, eta in signs/np.sqrt(3):
        sx, sy = signs.T
        N = (1+sx*xi)*(1+sy*eta)/4
        D = np.column_stack((sx*(1+sy*eta), sy*(1+sx*xi)))/4
        J = X.T@D
        Q.append((N, np.linalg.det(J), D@np.linalg.inv(J)))
    # Directional derivatives of psi from its explicit 2D quadratic formula.
    def dpsi(de):
        return (bulk*max(trace, 0.)*np.trace(de)
                + mu*((elastic[0, 0]-elastic[1, 1])*(de[0, 0]-de[1, 1])
                      +4*elastic[0, 1]*de[0, 1]))

    directions = {}
    perturb = np.array([[.02, .013], [.013, -.009]])
    vec = np.zeros_like(U)
    vec[:, [indices['ux'], indices['uy']]] = X@perturb.T
    directions['affine_strain'] = (vec.ravel(), dpsi(perturb), 0., 0.)
    for name in ['c', 'vfrac', 'd', 'q']:
        vec = np.zeros_like(U)
        vec[:, indices[name]] = 1.
        directions[name] = (vec.ravel(), dpsi(-p['beta']*np.eye(2)) if name == 'c' else 0.,
                            float(name == 'vfrac'), float(name == 'd'))
    results = {}
    cases = [('loading', .001, .001), ('unloading', .001, .1), ('fatigue_floor', .1, .001)]
    if relaxation:
        cases.append(('capped_relaxation', .001, .001))
    for label, alpha_old, W_old in cases:
        old = np.zeros((4, 4))
        old[:, 0] = [.0001, .0002, .0003, 1e16 if label == 'fatigue_floor' else 10.]
        if label == 'capped_relaxation':
            old[3, 0] = max(10., 2*cap)
        old[:, 1] = alpha_old
        old[:, 2] = W_old
        old[:, 3] = psi*np.array([.2, theta/2, 1.2, .3])
        rhs, K, state, pnewdt = direct(module, contract, X, U.ravel(),
                                       previous=previous.ravel(), state=old.ravel(),
                                       props=props, dt=args.dt)
        assert pnewdt == 1.
        expected_R = np.zeros(4)
        expected_Kv = {name: np.zeros(4) for name in directions}
        expected_state = []
        phase_history = []
        gamma_diagonal = np.zeros(4)
        gamma_matrix = np.zeros_like(K)
        gamma_old_damage_tangent = np.zeros_like(K)
        branches = []
        rows = np.arange(4)*dpn+indices['d']
        delta_damage = damage-previous[:, indices['d']]
        for gp, (N, weight, B) in enumerate(Q):
            H0, A0, W0, psi0 = old[gp]
            loading = W > W0
            alpha = A0+max(W-W0, 0.)
            if alpha <= p['alpha_T']:
                f, fp = 1., 0.
            else:
                z = max(0., 1-p['k_fg']*np.log10(alpha/p['alpha_T']))
                f = max(1e-15, z*z)
                fp = -2*p['k_fg']*z/(alpha*np.log(10.)) if f > 1e-15 else 0.
            phasepsi = min(psi, max(psi0, theta*psi))
            if psi < max(psi0, theta*psi):
                slope, branch = 1., 'current_cap'
            elif theta*psi > psi0:
                slope, branch = theta, 'predictor_floor'
            else:
                slope, branch = 0., 'accepted_energy'
            candidate = C*phasepsi/f
            active = candidate >= H0
            phaseH = max(H0, candidate)
            if args.relaxation_form == 'accepted_h':
                gamma = effective_relaxation*min(H0, cap)
                gamma_H = 0.
            elif args.relaxation_form == 'accepted_reaction':
                gamma = effective_relaxation*(1+min(a*H0, cap))
                gamma_H = 0.
            else:
                gamma = effective_relaxation*(1+min(a*phaseH, cap))
                gamma_H = effective_relaxation*a if a*phaseH < cap else 0.
            gamma_diagonal += N*weight*gamma
            physicalH = max(H0, C*psi/f)
            expected_state.append([physicalH, alpha, W, psi])
            phase_history.append(phaseH)
            branches.append(branch if active else 'frozen_H')
            expected_R += N*weight*(damage-a*vfrac*phaseH+gamma*delta_damage)

            def phase_derivative(psip, vp):
                alphap = (g*psip+2*a*vfrac*psi*vp) if loading else 0.
                return (C*slope*psip/f-candidate*fp*alphap/f) if active else 0.

            for name, (_, psip, vp, dp) in directions.items():
                Hp = phase_derivative(psip, vp)
                expected_Kv[name] += N*weight*((1+gamma)*dp-a*vp*phaseH
                                              +(-a*vfrac+gamma_H*delta_damage)*Hp)
            # Closed-form per-DOF H derivative isolates the complete viscosity
            # tangent, including its current-H chain rather than only a diagonal.
            Hgradient = np.zeros(K.shape[1])
            for node in range(4):
                bx, by = B[node]
                Hgradient[node*dpn+indices['ux']] = phase_derivative(
                    dpsi(np.array([[bx, .5*by], [.5*by, 0.]])), 0.)
                Hgradient[node*dpn+indices['uy']] = phase_derivative(
                    dpsi(np.array([[0., .5*bx], [.5*bx, by]])), 0.)
                Hgradient[node*dpn+indices['c']] = phase_derivative(
                    N[node]*dpsi(-p['beta']*np.eye(2)), 0.)
                Hgradient[node*dpn+indices['vfrac']] = phase_derivative(0., N[node])
            gamma_matrix[rows, rows] += N*weight*gamma
            gamma_matrix[rows, :] += np.outer(N*weight*gamma_H*delta_damage, Hgradient)
            gamma_old_damage_tangent[rows, :] += np.outer(N*weight*gamma_H, Hgradient)
        R = -rhs[rows]
        residual_error = float(np.max(abs(R-expected_R)/np.maximum(1., abs(expected_R))))
        state_error = float(np.max(abs(state.reshape(4, 4)-expected_state)
                                  /np.maximum(1., np.abs(expected_state))))
        assert residual_error < 1e-12 and state_error < 1e-12
        tangent_errors = {}
        for name, (direction, _, _, _) in directions.items():
            actual = (K@direction)[rows]
            expected = expected_Kv[name]
            error = float(np.max(abs(actual-expected)/np.maximum(1., abs(expected))))
            assert error < 1e-11, (label, name, error, actual, expected)
            tangent_errors[name] = error
        assert branches == ['accepted_energy', 'predictor_floor', 'current_cap', 'frozen_H']
        results[label] = dict(branches=branches, phase_residual_scaled_error=residual_error,
                              phase_history=phase_history,
                              relaxation_diagonal=gamma_diagonal.tolist(),
                              physical_history=np.asarray(expected_state)[:, 0].tolist(),
                              history_scaled_error=state_error, tangent_scaled_errors=tangent_errors)
        if label == 'loading':
            # Removing the old-energy state must change the phase equation,
            # while committed physical history remains independent of it.
            erased = old.copy()
            erased[:, 3] = 0.
            bad = direct(module, contract, X, U.ravel(), previous=previous.ravel(),
                         state=erased.ravel(), props=props, dt=args.dt)
            difference = float(np.max(abs(bad[0][rows]-rhs[rows])))
            assert difference > .01 and np.array_equal(bad[2], state)
            results['erased_accepted_psi_control'] = dict(
                phase_residual_difference=difference,
                physical_history_bitwise_unchanged=True)
            # Calls discarded by the host must leave this state fully repeatable.
            for factor in (.3, 2., .7):
                trial = U.copy()
                trial[:, [indices['ux'], indices['uy']]] *= factor
                direct(module, contract, X, trial.ravel(), previous=previous.ravel(),
                       state=old.ravel(), props=props, dt=args.dt)
            repeat = direct(module, contract, X, U.ravel(), previous=previous.ravel(),
                            state=old.ravel(), props=props, dt=args.dt)
            assert all(np.array_equal(x, y) for x, y in zip((rhs, K, state, pnewdt), repeat))
            results['discarded_trials'] = 'bitwise identical RHS, K, SVARS and PNEWDT'
        if relaxation and label in ('loading', 'capped_relaxation'):
            # Independently isolate the entire postgenerated viscosity term.
            zero_props = props.copy()
            zero_props[contract['prop_names'].index('phase_relaxation_factor')] = 0.
            zero = direct(module, contract, X, U.ravel(), previous=previous.ravel(),
                          state=old.ravel(), props=zero_props, dt=args.dt)
            matrix_error = float(np.max(abs(K-zero[1]-gamma_matrix)))
            residual_added = gamma_diagonal*(damage-previous[:, indices['d']])
            residual_added_error = float(np.max(abs(-rhs[rows]+zero[0][rows]-residual_added)))
            assert matrix_error < 1e-11 and residual_added_error < 1e-11
            assert np.array_equal(state, zero[2])
            # Change only U-DU, leaving physical fields and all old SVARS fixed.
            changed = previous.copy()
            delta_old = .03*np.array([1., -2., 3., -.5])
            changed[:, indices['d']] += delta_old
            shifted = direct(module, contract, X, U.ravel(), previous=changed.ravel(),
                             state=old.ravel(), props=props, dt=args.dt)
            old_damage_error = float(np.max(abs(-shifted[0][rows]+rhs[rows]
                                               +gamma_diagonal*delta_old)))
            predicted_tangent_shift = np.zeros_like(K)
            predicted_tangent_shift[rows, :] = -delta_old[:, None]*gamma_old_damage_tangent[rows, :]
            old_damage_tangent_error = float(np.max(abs(shifted[1]-K-predicted_tangent_shift)))
            assert old_damage_error < 1e-11
            assert old_damage_tangent_error < 1e-11 and np.array_equal(state, shifted[2])
            results[label]['isolated_relaxation'] = dict(
                tangent_matrix_absolute_error=matrix_error,
                residual_absolute_error=residual_added_error,
                accepted_old_damage_absolute_error=old_damage_error,
                accepted_old_damage_tangent_absolute_error=old_damage_tangent_error,
                physical_history_bitwise_unchanged=True,
                old_damage_shift_tangent_bitwise_unchanged=bool(np.array_equal(K, shifted[1])))
    source = contract['source']+['element/uel_contract.json']
    results['source_sha256'] = {name: hashlib.sha256((root/name).read_bytes()).hexdigest() for name in source}
    results['theta'] = theta
    results['relaxation_factor'] = relaxation
    results['effective_relaxation_factor'] = effective_relaxation
    results['relaxation_cap'] = cap
    results['relaxation_form'] = args.relaxation_form
    results['dt'] = args.dt
    results['scope'] = 'independent closed-form uniform-field phase row, mixed per-GP history branches'
    args.out.write_text(json.dumps(results, indent=2)+'\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
