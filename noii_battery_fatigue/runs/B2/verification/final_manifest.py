"""Gate the final evidence record on completed, source-matching host studies."""
from pathlib import Path
import json
import hashlib
import subprocess

ROOT=Path(__file__).resolve().parents[1]


if __name__=='__main__':
    def read(name):return json.loads((ROOT/'verification'/name).read_text())
    r=read('results.json');audit=read('output_audit.json');provenance=read('provenance.json')
    assert not r['pending_final_source_runs'],r['pending_final_source_runs']
    assert not r['failed_final_source_runs'],r['failed_final_source_runs']
    assert all(v=='passed' for v in r['tests'].values()),r['tests']
    assert audit['status']=='passed'
    runs=r['completed_final_source_runs']
    assert len(runs)>=21,'The planned final-source study cohort is incomplete'
    assert len(r['tests'])>=11,'A validation gate is missing'
    assert {x['path'] for x in runs}=={x['path'] for x in audit['runs']},'Audit does not cover all completed runs'
    sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in r['source_sha256']}
    assert sources==r['source_sha256']==provenance['source_sha256']
    contract_sha=hashlib.sha256((ROOT/'element/uel_contract.json').read_bytes()).hexdigest()
    assert contract_sha==provenance['contract_sha256']
    for name,expected in provenance['python_source_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected,name
    for item in runs:
        assert item['accepted_output_steps']==item['cycles']*item['steps_per_cycle']
        assert all(v>=0 for v in item['history_min_change'].values())
    for name in ('element/uel_contract.json','design.md','research_log.md','STATUS.md','verification/README.md','verification/self_review.md'):
        assert (ROOT/name).is_file() and (ROOT/name).stat().st_size>0,name
    for rec in r['field_comparisons']:
        assert rec['interpolation']['status']=='passed'
        assert rec['comparison']['comparison_version']==2
    from update_evidence import plans
    expected_comparisons={(prod,ref,cycle) for _,prod,ref,cycle in plans()}
    actual_comparisons={(rec['comparison']['production'],rec['comparison']['reference'],
                         rec['comparison']['cycle']) for rec in r['field_comparisons']}
    assert expected_comparisons<=actual_comparisons,'A planned field comparison is missing'
    refs=r['independent_references']
    assert all(x.get('output_coverage',{}).get('status')=='passed' for x in refs if x['status']=='ok')
    failed_refs=[x['path'] for x in refs if x['status']!='ok']
    result=dict(status='passed',scope='Completeness and integrity of submitted-element evidence; reference limits are listed separately',
        recorded_at=subprocess.check_output(['date','-u','+%Y-%m-%dT%H:%M:%SZ'],text=True).strip(),
        completed_global_runs=len(runs),accepted_nominal_outputs=sum(x['accepted_output_steps'] for x in runs),
        rejected_global_attempts=sum(x.get('rejected_attempts',0) for x in runs),
        audited_scalar_values=sum(x['checked_scalar_values'] for x in audit['runs']),
        passed_check_records=r['tests'],field_comparisons=len(r['field_comparisons']),
        planned_field_comparisons_complete=len(expected_comparisons),
        completed_independent_references=sum(x['status']=='ok' for x in refs),
        failed_independent_references=failed_refs,
        pending_independent_references=r.get('pending_independent_references',[]),
        source_sha256=sources,contract_sha256=contract_sha,
        required_deliverables_present=True,
        known_limits=['Raw projected pressure gradient does not converge strongly in the distorted manufactured family',
            'Approximately first-order manufactured concentration convergence',
            'First-order accepted-history time error','Coarse Q1 host energy reconstruction error',
            'Consistent-mass concentration can undershoot in unresolved low-mobility layers',
            'Distorted-mesh damage bounds have only an empirical validation envelope'],
        unverified=['Hidden problems and organizer reference','Abaqus execution',
            'Mesh/time resolutions beyond the measured envelope'])
    (ROOT/'verification/final_manifest.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
