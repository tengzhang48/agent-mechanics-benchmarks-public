"""Independent small-matrix oracle for extreme heterogeneous frozen H.

The exact nodal lumped-reaction system is (ell^2 L+diag(m+a*S))*v=m.
The unit-square Q1 stiffness is entered analytically, not obtained from the
UEL or generated code. The host solves d,q,v from its zero auxiliary guess.
"""
import json

import numpy as np

from common import ROOT, compiled, fehost, source_identity


def main():
    module, contract = compiled()
    params = dict(zip(contract['prop_names'], contract['props']))
    a, ell = 1-params['kappa'], params['ell']
    local_L = np.array([[4., -1., -2., -1.], [-1., 4., -1., -2.],
                        [-2., -1., 4., -1.], [-1., -2., -1., 4.]])/6
    signs = np.array([[-1., -1.], [1., -1.], [1., 1.], [-1., 1.]])
    N = np.array([(1+signs[:, 0]*xi)*(1+signs[:, 1]*eta)/4
                  for xi, eta in signs/np.sqrt(3)])
    records = []
    cases = [
        ('single_high_gp', [[1e15, 1., 1., 1.]]),
        ('alternating_high_gp', [[1e15, 1., 1e15, 1.]]),
        ('uniform_high', [[1e15]*4]),
        ('two_element_jump', [[1e15]*4, [1.]*4]),
    ]
    for name, histories in cases:
        histories = np.asarray(histories)
        ne = len(histories)
        if ne == 1:
            nodes = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
            elems = np.array([[0, 1, 2, 3]])
        else:
            nodes = np.array([[0., 0.], [1., 0.], [2., 0.], [0., 1.], [1., 1.], [2., 1.]])
            elems = np.array([[0, 1, 4, 3], [1, 2, 5, 4]])
        m, S = np.zeros(len(nodes)), np.zeros(len(nodes))
        L = np.zeros((len(nodes), len(nodes)))
        for ids, H in zip(elems, histories):
            L[np.ix_(ids, ids)] += local_L
            m[ids] += .25
            S[ids] += N.T@H/4
        reference_v = np.linalg.solve(ell*ell*L+np.diag(m+a*S), m)
        assert np.all(reference_v > 0.)
        element = fehost.ContractUEL(module, contract, ne, .2)
        hcols = np.asarray(contract['history_slots']['H'])-1
        element.svars[:, hcols] = histories
        dpn = element.dof_per_node
        initial = np.zeros((len(nodes), dpn))
        initial[:, element.comp('c')] = params['c0']
        bc = {dpn*i+element.comp(field): params['c0'] if field == 'c' else 0.
              for i in range(len(nodes)) for field in ['ux', 'uy', 'c']}
        solution, report = fehost.solve(nodes, elems, element, lambda _: bc, 1, .2,
                         U0=initial.ravel(), tol=1e-10,
                         residual_scale=fehost.residual_weights(element, contract))
        values = solution.reshape(len(nodes), dpn)
        damage, intact = values[:, element.comp('d')], values[:, element.comp('vfrac')]
        err = float(np.max(abs(intact-reference_v)/reference_v))
        constraint = float(np.max(abs((damage-1)+intact)))
        assert err < 1e-10 and constraint < 1e-12 and report.cutbacks == 0
        assert np.array_equal(element.svars[:, hcols], histories)
        records.append(dict(case=name, H=histories.tolist(),
                            v_reference=reference_v.tolist(), v_observed=intact.tolist(),
                            v_max_relative_error=err, intact_constraint_max=constraint,
                            d_min=float(damage.min()), d_max=float(damage.max()),
                            H_bitwise_preserved=True, solver=report.summary()))
    result = dict(scope='Independent analytical unit-cell matrix; frozen heterogeneous H, not a fatigue path',
                  source_sha256=source_identity(), results=records)
    (ROOT/'verification/heterogeneous_fixed_history.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
