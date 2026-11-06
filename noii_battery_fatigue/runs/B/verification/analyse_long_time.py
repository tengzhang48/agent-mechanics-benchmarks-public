"""Source-matched temporal refinement through 20 D0 cycles with active fatigue."""
import argparse
import hashlib
import json

import numpy as np

from common import ROOT, source_identity
from analyse_studies import compare_fields, independent_observables, write_csv


def polyline_integrals(times, values):
    """Return exact integral(abs(y)) and integral(y**2) for linear pieces."""
    h = np.diff(times)
    left, right = values[:-1], values[1:]
    magnitude = abs(left) + abs(right)
    area = h*magnitude/2
    crosses = (left < 0) & (right > 0) | (left > 0) & (right < 0)
    area[crosses] = h[crosses]*(left[crosses]**2+right[crosses]**2)/(2*magnitude[crosses])
    square = h*(left**2+left*right+right**2)/3
    return float(area.sum()), float(square.sum())


def time_norms(coarse, fine, c0):
    """Integrate piecewise-linear nominal-step curves, including known t=0."""
    initial=dict(d_max=0.,D_int=0.,C_int=c0,E_el=0.,A_damage=0.)
    times=np.unique(np.r_[0.,coarse['t'],fine['t']])
    out={}
    for key,start in initial.items():
        reference=np.interp(times,np.r_[0.,fine['t']],np.r_[start,fine[key]])
        values=np.interp(times,np.r_[0.,coarse['t']],np.r_[start,coarse[key]])
        difference=values-reference
        norm1,norm2=polyline_integrals(times,reference)
        error1,error2=polyline_integrals(times,difference)
        out[key+'_relative_time_l1']=error1/norm1 if norm1 else None
        out[key+'_relative_time_l2']=float(np.sqrt(error2/norm2)) if norm2 else None
        maximum=float(np.max(abs(difference)))
        reference_peak=float(np.max(abs(reference)))
        out[key+'_max_abs_diff']=maximum
        out[key+'_max_rel_diff']=maximum/reference_peak if reference_peak else None
    crossing=np.flatnonzero(coarse['d_max']>=.9)
    if len(crossing):
        i=int(crossing[0])
        before_t,before_d=(coarse['t'][i-1],coarse['d_max'][i-1]) if i else (0.,0.)
        out['dmax_09_first_crossing_linear_time']=float(before_t+(.9-before_d)*
                   (coarse['t'][i]-before_t)/(coarse['d_max'][i]-before_d))
    else:
        out['dmax_09_first_crossing_linear_time']=None
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--through-spc',type=int,choices=[1280,2560],default=2560)
    parser.add_argument('--interim',action='store_true')
    args=parser.parse_args()
    assert args.interim or args.through_spc==2560, 'Only the full cohort writes canonical results.'
    identity = source_identity()
    contract = json.loads((ROOT/'element/uel_contract.json').read_text())
    params = json.loads((ROOT/'host/dev_problem/d0.json').read_text())['parameters']
    Hcols = np.asarray(contract['history_slots']['H'])-1
    Acols = np.asarray(contract['history_slots']['alpha'])-1
    cases = []
    audits = []
    artifact_records = []
    for spc in [value for value in [40,80,160,320,640,1280,2560] if value<=args.through_spc]:
        version = 'selected' if spc == 40 else 'long_time'
        label = f'D0_structured_16_s1_c20_spc{spc}'
        folder = ROOT/f'verification/runs/{version}/{label}'
        run = json.loads((folder/'run.json').read_text())
        assert run['status'] == 'ok' and run['cycles'] == 20
        assert run['accepted_output_steps'] == 20*spc
        for name, expected in run['element_source_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected
        assert all(run['element_source_sha256'][name] == expected for name, expected in identity.items())
        observables = np.atleast_1d(np.genfromtxt(folder/'observables.csv', delimiter=',', names=True))
        assert len(observables) == 20*spc
        assert np.all(np.diff(observables['t']) > 0)
        assert np.isfinite(observables.tolist()).all()
        assert all(value >= 0 for value in run['history_min_change'].values())
        expected_names={f'fields_cycle_{cycle:03d}.npz' for cycle in range(1,21)}
        assert {file.name for file in folder.glob('fields_cycle_*.npz')} == expected_names
        for cycle in range(1,21):
            fields = dict(np.load(folder/f'fields_cycle_{cycle:03d}.npz'))
            assert fields['svars'].shape == (256,16) and fields['U'].size == 289*6
            assert np.isfinite(fields['U']).all() and np.isfinite(fields['svars']).all()
            computed = independent_observables(fields, params)
            stored = observables[cycle*spc-1]
            bridge = max(abs(computed[key]-stored[key])/max(abs(computed[key]),1e-12) for key in computed)
            assert bridge < 2e-9
            audits.append(dict(label=label,spc=spc,cycle=cycle,bridge_relative_error=bridge))
        artifact_records.append(dict(path=str(folder.relative_to(ROOT)),
            artifact_sha256={file.name:hashlib.sha256(file.read_bytes()).hexdigest()
                             for file in sorted(folder.iterdir()) if file.is_file()
                             and file.suffix in ('.csv','.json','.npz')}))
        cases.append((spc, label, folder, run, fields, observables))
    ref = cases[-1]
    rows = []
    for spc, label, folder, run, fields, observables in cases:
        alpha = float(fields['svars'][:, Acols].max())
        f = 1. if alpha <= params['alpha_T'] else max(1e-15,max(0.,1-params['k_fg']*np.log10(alpha/params['alpha_T']))**2)
        row = dict(spc=spc, label=label, path=str(folder.relative_to(ROOT)), reference=ref[1],
                   cutbacks=run['cutbacks'], H_max_final=float(fields['svars'][:, Hcols].max()),
                   alpha_max_final=alpha, f_min_final=float(f),
                   **time_norms(observables,ref[-1],params['c0']),
                   **compare_fields(fields, ref[-2]))
        for name, cols in [('H',Hcols),('alpha',Acols),('W',Hcols+2),('psi_prev',Hcols+3)]:
            exact = ref[-2]['svars'][:,cols]
            row[name+'_final_relative_l2'] = float(np.linalg.norm(fields['svars'][:,cols]-exact)/max(np.linalg.norm(exact),1e-30))
        rows.append(row)
    successive = [dict(coarse_spc=coarse[0], fine_spc=fine[0],
                       **time_norms(coarse[-1],fine[-1],params['c0']),
                       **compare_fields(coarse[-2],fine[-2]))
                  for coarse,fine in zip(cases[:-1],cases[1:])]
    result = dict(scope='Same-mesh 20-cycle temporal self-convergence with active fatigue; no independent full coupled reference',
                  source_sha256=identity, reference_spc=args.through_spc, records=rows,
                  final_planned_cohort_complete=args.through_spc==2560,
                  successive_refinements=successive,
                  field_snapshots_audited=len(audits),
                  curve_norm_definition='Exact interval integrals and maximum differences of piecewise-linear nominal output curves, including exact D0 initial values, on the union of both knot sets. Absolute-value integrals split at sign changes. First dmax=.9 crossing is linearly interpolated, not an exact physical jump time.',
                  max_bridge_relative_error=max(row['bridge_relative_error'] for row in audits),
                  case_artifacts=artifact_records)
    out = ROOT/'verification/long_time_refinement'
    out.mkdir(exist_ok=True)
    suffix=f'_interim_to{args.through_spc}' if args.interim else ''
    write_csv(out/f'comparison{suffix}.csv', rows)
    write_csv(out/f'successive_comparison{suffix}.csv', successive)
    write_csv(out/f'field_audit{suffix}.csv', audits)
    result_name=f'interim_to{args.through_spc}.json' if args.interim else 'results.json'
    (out/result_name).write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
