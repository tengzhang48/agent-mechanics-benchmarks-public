"""Global fixed-H conditioning and exact relaxation recurrence on the host."""
import json

import numpy as np

from common import ROOT, compiled, fehost, source_identity


def main():
    module, contract = compiled()
    X = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
    elems = np.array([[0, 1, 2, 3]], int)
    dpn = contract['dof_per_node']
    params = dict(zip(contract['prop_names'], contract['props']))
    eta = params.get('phase_relaxation_factor', 0.)
    cap = params.get('phase_relaxation_cap', 1000.)
    a = 1-params['kappa']
    checks = []
    for H in [1., 1e8, 1e12, 1e15, 1e20]:
        el = fehost.ContractUEL(module, contract, 1, .2)
        slots = np.asarray(contract['history_slots']['H'])-1
        el.svars[:, slots] = H
        initial = np.zeros(4*dpn)
        initial[el.comp('c')::dpn] = params['c0']
        gamma = eta*min(H, cap)
        previous_d = 0.
        errors = []
        sequence = []

        def bc(fraction):
            return {dpn*i+el.comp(name): (params['c0'] if name == 'c' else 0.)
                    for i in range(4) for name in ['ux', 'uy', 'c']}

        def observe(step, fraction, U, metrics=None):
            nonlocal previous_d
            # Independent scalar recurrence using the old physical d; this
            # remains meaningful when 1-d has rounded to zero at huge H.
            v_expected = (1+gamma*(1-previous_d))/(1+gamma+a*H)
            d_expected = 1-v_expected
            d, v = U[el.comp('d')::dpn], U[el.comp('vfrac')::dpn]
            error_d = float(np.max(abs(d-d_expected)))
            error_v = float(np.max(abs(v-v_expected))/max(v_expected, 1e-300))
            assert error_d < 2e-10 and error_v < 2e-7, (H, step, error_d, error_v)
            errors.append((error_d, error_v))
            sequence.append(dict(step=step, d=float(d.mean()), v=float(v.mean()),
                                 exact_d=d_expected, exact_v=v_expected))
            previous_d = float(d.mean())

        U, report = fehost.solve(X, elems, el, bc, 24, .2, U0=initial,
                                tol=1e-10, step_callback=observe,
                                residual_scale=fehost.residual_weights(el, contract))
        target_v = 1/(1+a*H)
        assert np.max(abs(U[el.comp('d')::dpn]-(1-target_v))) < 2e-10
        checks.append(dict(H=H, gamma=gamma, status='ok',
                           steady_d=1-target_v, steady_v=target_v,
                           final_d_min=float(U[el.comp('d')::dpn].min()),
                           final_d_max=float(U[el.comp('d')::dpn].max()),
                           final_v_min=float(U[el.comp('vfrac')::dpn].min()),
                           max_recurrence_d_error=max(x[0] for x in errors),
                           max_recurrence_v_relative_error=max(x[1] for x in errors),
                           recurrence=sequence, **report.summary()))
    result = dict(element_source_sha256=source_identity(), checks=checks)
    (ROOT/'verification/fixed_history_global.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
