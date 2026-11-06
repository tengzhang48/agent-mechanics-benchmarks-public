"""Audit the delivered source identities and the selected evidence cohort.

This validates artifact completeness and provenance, not hidden-case accuracy.
Run analyse_studies.py first. --allow-pending is for an explicit interim audit.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

from common import ROOT, source_identity


def read(name):
    return json.loads((ROOT/name).read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def public_plan():
    """The fixed selected cohort; no solver or runner is invoked."""
    cases = [(mesh, n, 5, 80) for n in (8, 16, 32, 64)
             for mesh in ('structured', 'distorted')]
    cases += [('structured', 16, 5, spc) for spc in (10, 20, 40, 160, 320, 640, 1280)]
    cases += [('structured', 32, 5, spc) for spc in (40, 160)]
    cases += [(mesh, 16, 20, 40) for mesh in ('structured', 'distorted')]
    cases += [('structured', 32, 20, 80)]
    plan = {f'D0_{mesh}_{n}_s1_c{cycles}_spc{spc}':
            dict(mesh=mesh, n=n, cycles=cycles, steps_per_cycle=spc, seed=1)
            for mesh, n, cycles, spc in cases}
    assert len(plan) == 20 and sum(case['cycles'] for case in plan.values()) == 145
    return plan


def check_hashes(folder, manifest, context):
    for file, expected in manifest.items():
        path = (folder/file).resolve()
        assert path.is_relative_to(ROOT), ('outside-workspace manifest', context, file)
        assert sha(path) == expected, ('changed recorded artifact', context, file)


def artifact_hashes(folder):
    """Record the durable scientific files, excluding compiler/cache output."""
    return {file.name: sha(file) for file in sorted(folder.iterdir())
            if file.is_file() and (file.suffix in ('.csv', '.json', '.npz')
                                   or file.name == 'stdout.txt')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--allow-pending', action='store_true')
    args = parser.parse_args()
    identity = source_identity()
    contract = read('element/uel_contract.json')
    compiled_sources = contract['source']
    assert contract['dof_names'] == ['ux', 'uy', 'c', 'd', 'q', 'vfrac']
    assert contract['nsvars'] == 16
    assert contract['history_slots'] == {'H': [1, 5, 9, 13], 'alpha': [2, 6, 10, 14]}
    assert len(contract['props']) == len(contract['prop_names']) == 12
    for name in ['design.md', 'STATUS.md', 'research_log.md', 'verification/README.md']:
        assert (ROOT/name).is_file()
    provenance = read('verification/provenance.json')
    assert provenance['unchanged_host_and_tool_audit']['host_identical']
    assert provenance['unchanged_host_and_tool_audit']['generator_identical']
    check_hashes(ROOT, provenance['host_sha256'], 'live host source')
    check_hashes(ROOT, provenance['generator_sha256'], 'live generator source')
    for name, expected in provenance['element_sha256'].items():
        assert sha(ROOT/name) == expected, ('changed element since provenance', name)
    assert read('verification/generation.json')['sha256'] == identity[compiled_sources[0]]

    exact_files = [
        'verification/element_checks.json', 'verification/transactions.json',
        'verification/fixed_history_global.json', 'verification/pressure_gradient_patch.json',
        'verification/heterogeneous_fixed_history.json',
        'verification/uniform_history_oracle.json', 'verification/psilag_phase_oracle.json',
        'verification/degradation_gradient.json',
        'verification/constrained_oracle_results/results.json',
        'verification/quadratic_drift_results.json',
        'verification/manufactured_history_results/results.json',
    ]
    exact_records = []
    for name in exact_files:
        data = read(name)
        manifest = data.get('element_source_sha256', data.get('source_sha256', {}))
        assert set(identity) <= set(manifest), ('missing source/contract identity', name)
        for file, expected in manifest.items():
            assert sha(ROOT/file) == expected, ('stale source evidence', name, file)
        exact_records.append(dict(path=name, sha256=sha(ROOT/name),
                                  contract_identity_recorded='element/uel_contract.json' in manifest))
    transactions = read('verification/transactions.json')
    for key in ['bitwise_identical_final_U', 'bitwise_identical_final_svars',
                'bitwise_identical_all_output_fields', 'bitwise_identical_accepted_times']:
        assert transactions[key] is True
    direct = read('verification/element_checks.json')
    assert direct['zero_time']['finite'] and direct['zero_time']['state_preserved']
    assert direct['zero_time']['pnewdt'] < 1

    pressure_path = ROOT/'verification/pressure_floor_diagnostic/results.json'
    pressure = read(str(pressure_path.relative_to(ROOT)))
    assert pressure['source_sha256'] == identity
    check_hashes(ROOT, pressure['input_sha256'], 'pressure diagnostic inputs')
    # This diagnostic records workspace-relative artifact paths; resolve them
    # against ROOT, then enforce and hash within its own artifact directory.
    pressure_folder = pressure_path.parent.resolve()
    pressure_artifacts = {}
    for file, expected in pressure['artifact_sha256'].items():
        path = (ROOT/file).resolve()
        assert path.is_relative_to(pressure_folder), ('outside pressure diagnostic', file)
        pressure_artifacts[str(path.relative_to(pressure_folder))] = expected
    check_hashes(pressure_folder, pressure_artifacts, 'pressure diagnostic artifacts')

    selected = read('verification/custom_runs/selected/source_bound_cases.json')
    for file in identity:
        assert selected['source_sha256'][file] == identity[file]
    assert not selected['pending']
    assert len(selected['cases']) == len({case['case'] for case in selected['cases']}) == 15
    check_hashes(ROOT/selected['source_snapshot'], selected['source_sha256'], 'authored source snapshot')
    selected_records = []
    failures = []
    artifacts_checked = 0
    for case in selected['cases']:
        folder = ROOT/case['durable_path']
        data = read(str(folder.relative_to(ROOT)/'run.json'))
        record = read(str(folder.relative_to(ROOT)/'provenance.json'))
        assert data['status'] == case['status'] == record['summary']['status']
        for file in identity:
            assert record['source_sha256'][file] == identity[file]
        check_hashes(folder, record['artifact_sha256'], case['case'])
        artifacts_checked += len(record['artifact_sha256'])
        selected_records.append(dict(case=case['case'], status=data['status'],
                                     path=case['durable_path']))
        if data['status'] != 'ok':
            failures.append(case['case'])
            assert (folder/'failed_trial.npz').is_file()
    assert failures == ['spatial_floor_distorted16_s3_c2_spc80'], failures

    alternatives = []
    rejected_relaxations = []
    alternative_artifacts_checked = 0
    for family, expected_failures in [('accepted_history_viscosity', 3),
                                      ('accepted_stiffness_viscosity', 4),
                                      ('current_stiffness_viscosity', 3),
                                      ('sqrt_dt_stiffness_viscosity', 3)]:
        folder = ROOT/'verification/alternatives'/family
        outcomes = read(str((folder/'outcomes.json').relative_to(ROOT)))
        frozen = read(str((folder/'verification/frozen_source.json').relative_to(ROOT)))
        # Some generic verification/runner files were fingerprinted but not
        # archived. All actual element source inputs must be present and match.
        source_files = {file: data['sha256'] for file, data in frozen['source_files'].items()
                        if file.startswith('element/')}
        assert set(identity) <= set(source_files)
        check_hashes(folder, source_files, family+' source snapshot')
        assert len(outcomes) == len({record['path'] for record in outcomes})
        family_failures = 0
        for record in outcomes:
            run_path = folder/record['path']
            data = read(str(run_path.relative_to(ROOT)))
            assert data['status'] == record['status'] == record['summary']['status']
            check_hashes(run_path.parent, record['artifact_sha256'], family+'/'+record['path'])
            alternative_artifacts_checked += len(record['artifact_sha256'])
            assert len(record['runtime_prop_names']) == len(record['runtime_props']) == 14
            assert record['runtime_props_evidence']
            alternatives.append(dict(family=family, path=str(run_path.relative_to(ROOT)),
                                     status=data['status']))
            if data['status'] != 'ok':
                assert (run_path.parent/'failed_trial.npz').is_file()
                family_failures += 1
                rejected_relaxations.append(str(run_path.parent.relative_to(ROOT)))
        assert family_failures == expected_failures, (family, family_failures)
    assert len(alternatives) == 26 and len(rejected_relaxations) == 13

    quad_root = ROOT/'verification/quadrature_comparison/selected/source_snapshot'
    quad_provenance = read(str(quad_root.relative_to(ROOT)/'provenance.json'))
    quad_identity = quad_provenance['element_source_sha256']
    for file, expected in quad_identity.items():
        assert sha(quad_root/file) == expected
    assert quad_provenance['n_gauss_points'] == 9
    quad_contract = json.loads((quad_root/'element/uel_contract.json').read_text())
    assert quad_contract['nsvars'] == 36
    assert quad_contract['dof_names'] == contract['dof_names']
    assert sha(ROOT/'element/model.py') == quad_provenance['base_model_sha256']
    assert sha(ROOT/'element/lump_reaction.py') == quad_provenance['base_extension_sha256']
    quad_records = []
    for n in [16, 32]:
        name = f'verification/quadrature_comparison/selected/n{n}_3x3/run.json'
        data = read(name)
        assert data['status'] == 'ok' and data['element_source_sha256'] == quad_identity
        assert data['n'] == n and data['cycles'] == 3
        expected_fields = {f'fields_cycle_{cycle:03d}.npz' for cycle in range(1, 4)}
        folder = ROOT/name
        assert {file.name for file in folder.parent.glob('fields_cycle_*.npz')} == expected_fields
        quad_records.append(dict(n=n, status=data['status'], cutbacks=data['cutbacks'],
                                 artifact_sha256=artifact_hashes(folder.parent)))

    analysis = read('verification/study_analysis.json')
    assert analysis['version'] == 'selected' and analysis['element_source_sha256'] == identity
    plan = public_plan()
    public = []
    observed_snapshots = set()
    newest_field_mtime = 0.
    for path in sorted((ROOT/'verification/runs/selected').glob('*/run.json')):
        data = json.loads(path.read_text())
        assert path.parent.name in plan, ('unplanned selected public case', path.parent.name)
        for key, expected in plan[path.parent.name].items():
            assert data[key] == expected, ('case metadata mismatch', path.parent.name, key)
        assert set(identity) <= set(data['element_source_sha256'])
        for file, expected in data['element_source_sha256'].items():
            assert sha(ROOT/file) == expected, ('changed public source', path.parent.name, file)
        assert data['status'] == 'ok'
        expected_fields = {f'fields_cycle_{cycle:03d}.npz' for cycle in range(1, data['cycles']+1)}
        actual_fields = {file.name for file in path.parent.glob('fields_cycle_*.npz')}
        assert actual_fields == expected_fields, ('incomplete public field copy', path.parent.name)
        assert data['accepted_output_steps'] == data['cycles']*data['steps_per_cycle']
        with (path.parent/'observables.csv').open() as stream:
            assert sum(1 for _ in csv.DictReader(stream)) == data['accepted_output_steps']
        for cycle in range(1, data['cycles']+1):
            observed_snapshots.add((path.parent.name, cycle))
            newest_field_mtime = max(newest_field_mtime,
                                    (path.parent/f'fields_cycle_{cycle:03d}.npz').stat().st_mtime)
        public.append(dict(case=path.parent.name, cutbacks=data['cutbacks'],
                           run_json_sha256=sha(path), artifact_sha256=artifact_hashes(path.parent)))
    assert len(public) == analysis['completed_runs'] == analysis['ok_runs']
    assert len(observed_snapshots) == analysis['output_files_audited']
    bridge_records = []
    if observed_snapshots:
        with (ROOT/'verification/output_audit.csv').open() as stream:
            bridge_records = list(csv.DictReader(stream))
    assert len(bridge_records) == len(observed_snapshots)
    assert {(record['label'], int(record['cycle'])) for record in bridge_records} == observed_snapshots
    assert all(int(record['missing']) == int(record['nonfinite']) == 0 for record in bridge_records)
    assert all(float(record['bridge_max_relative_error']) < 2e-9 for record in bridge_records)
    assert (ROOT/'verification/study_analysis.json').stat().st_mtime >= newest_field_mtime, \
        'public snapshots changed since analysis; rerun analyse_studies.py'
    missing_public = sorted(set(plan)-{record['case'] for record in public})
    complete = not missing_public and len(observed_snapshots) == 145
    if not args.allow_pending:
        assert complete, 'selected public cohort incomplete'
    if complete:
        assert analysis['time_reference'] == 'D0_structured_16_s1_c5_spc1280'
    analysis_names = ['study_analysis.json', 'study_summary.csv', 'output_audit.csv',
                      'mesh_convergence.csv', 'field_convergence.csv',
                      'time_convergence.csv', 'time_field_convergence.csv',
                      'history_time_convergence.csv', 'quadrature_comparison.csv']
    analysis_artifacts = {name: sha(ROOT/'verification'/name) for name in analysis_names
                          if (ROOT/'verification'/name).is_file()}
    absent_analysis = sorted(set(analysis_names)-set(analysis_artifacts))
    if complete:
        assert not absent_analysis, ('missing completed-cohort analysis', absent_analysis)
    primary_complete=complete

    # Supplemental time refinement checks sustained fatigue at fixed n16;
    # these six cases are separate from the fixed 20-case primary cohort.
    long_plan={f'D0_structured_16_s1_c20_spc{spc}':spc for spc in (80,160,320,640,1280,2560)}
    long_records=[]
    for path in sorted((ROOT/'verification/runs/long_time').glob('*/run.json')):
        data=json.loads(path.read_text())
        assert path.parent.name in long_plan
        assert data['status']=='ok' and data['mesh']=='structured' and data['n']==16 and data['seed']==1
        assert data['cycles']==20 and data['steps_per_cycle']==long_plan[path.parent.name]
        assert data['accepted_output_steps']==20*data['steps_per_cycle']
        check_hashes(ROOT,data['element_source_sha256'],path.parent.name+' source')
        assert all(data['element_source_sha256'][file]==expected for file,expected in identity.items())
        assert {file.name for file in path.parent.glob('fields_cycle_*.npz')} == \
               {f'fields_cycle_{cycle:03d}.npz' for cycle in range(1,21)}
        long_records.append(dict(case=path.parent.name,cutbacks=data['cutbacks'],
                                 artifact_sha256=artifact_hashes(path.parent)))
    long_pending=sorted(set(long_plan)-{row['case'] for row in long_records})
    long_results_path=ROOT/'verification/long_time_refinement/results.json'
    reflection=read('verification/long_time_refinement/cycle16_reflection.json')
    assert reflection['source_sha256']==identity
    check_hashes(ROOT,reflection['input_sha256'],'cycle16 reflection diagnostic')
    long_comparison=None
    long_analysis_artifacts={}
    if not long_pending and long_results_path.exists():
        long_comparison=json.loads(long_results_path.read_text())
        assert long_comparison['source_sha256']==identity
        assert long_comparison['field_snapshots_audited']==140
        assert long_comparison['max_bridge_relative_error']<2e-9
        assert [row['spc'] for row in long_comparison['records']]==[40,80,160,320,640,1280,2560]
        assert long_comparison['reference_spc']==2560 and long_comparison['final_planned_cohort_complete']
        assert [(row['coarse_spc'],row['fine_spc']) for row in long_comparison['successive_refinements']]== \
               [(40,80),(80,160),(160,320),(320,640),(640,1280),(1280,2560)]
        for case in long_comparison['case_artifacts']:
            check_hashes(ROOT/case['path'],case['artifact_sha256'],'long-time comparison')
        expected_long_snapshots={(f'D0_structured_16_s1_c20_spc{spc}',cycle)
                                 for spc in (40,80,160,320,640,1280,2560) for cycle in range(1,21)}
        with (long_results_path.parent/'field_audit.csv').open() as stream:
            long_audits=list(csv.DictReader(stream))
        assert len(long_audits)==140
        assert {(row['label'],int(row['cycle'])) for row in long_audits}==expected_long_snapshots
        assert all(float(row['bridge_relative_error'])<2e-9 for row in long_audits)
        long_analysis_artifacts={name:sha(long_results_path.parent/name)
                                 for name in ('results.json','comparison.csv','successive_comparison.csv','field_audit.csv')}
    complete=primary_complete and not long_pending and long_comparison is not None
    if not args.allow_pending:
        assert complete,'supplemental long-cycle time-refinement evidence incomplete'
    result = dict(
        audited_utc=subprocess.check_output(['date', '-u', '+%Y-%m-%d %H:%M:%S UTC'], text=True).strip(),
        scope='Artifact/source completeness audit only; not organizer-condition certification or independent full coupled PDE validation',
        status='complete evidence cohort' if complete else 'explicitly interim',
        element_source_sha256=identity,
        unchanged_host_and_generator=True,
        perfect_organizer_condition_compliance_claimed=False,
        organizer_compliant_submission_valid=False,
        organizer_visibility_condition_breached=True,
        organizer_visibility_condition_reason=(
            'A broad ps listing exposed unrelated process command lines outside the workspace scope; '
            'ORGANIZER_NOTE.md states that such access invalidates the run.'),
        procedural_deviations=[
            'Two early public CLI builds created compiler scratch in verification/; removed and subsequent builds redirected to tmp/.',
            'A broad process-status listing exposed unrelated system command lines; no outside files or results used and the check was stopped.'
        ],
        exact_source_bound_records=exact_records,
        authored_cases=selected_records, authored_artifacts_rehashed=artifacts_checked,
        known_authored_failure=failures,
        rejected_relaxation_cases=rejected_relaxations,
        alternative_cases=alternatives,
        alternative_artifacts_rehashed=alternative_artifacts_checked,
        quadrature_cases=quad_records, public_cases=public,
        public_expected_runs=20, public_expected_field_snapshots=145,
        public_pending_cases=missing_public,
        primary_public_cohort_complete=primary_complete,
        supplemental_long_time_cases=long_records,
        supplemental_long_time_pending=long_pending,
        supplemental_comparison_sha256=sha(long_results_path) if long_comparison is not None else None,
        pressure_floor_diagnostic_results_sha256=sha(pressure_path),
        cycle16_reflection_sha256=sha(ROOT/'verification/long_time_refinement/cycle16_reflection.json'),
        supplemental_analysis_artifact_sha256=long_analysis_artifacts,
        public_analysis_artifact_sha256=analysis_artifacts,
        public_analysis_tables_not_yet_available=absent_analysis,
        public_field_snapshots_audited=analysis['output_files_audited'],
        max_public_bridge_relative_error=analysis['max_bridge_relative_error'],
        known_limitations=['heterogeneous fatigue-floor Newton failure',
                           'distorted-mesh recovered gradients can fail strong convergence',
                           'coarse-time notched fatigue activation is biased'],
        unverified=['Abaqus execution', 'organizer hidden cases',
                    'independent full coupled fracture reference',
                    'general nodal maximum principle',
                    'general heterogeneous fatigue-floor nonlinear robustness'])
    (ROOT/'verification/delivery_audit.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['audited_utc', 'status',
                    'authored_artifacts_rehashed', 'known_authored_failure',
                    'public_field_snapshots_audited']}, indent=2))


if __name__ == '__main__':
    main()
