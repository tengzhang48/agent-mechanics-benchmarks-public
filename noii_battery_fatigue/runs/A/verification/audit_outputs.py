"""Audit delivered artifacts, source provenance, finite fields and history."""
import hashlib
import json
from pathlib import Path
import numpy as np
from common import ROOT,write_json


def main():
    contract=json.loads((ROOT/'element'/'uel_contract.json').read_text())
    digest=hashlib.sha256((ROOT/contract['source']).read_bytes()).hexdigest()
    assert contract['dof_per_node']==len(contract['dof_names'])==len(contract['abaqus_dofs'])==6
    assert len(contract['props'])==len(contract['prop_names'])==11
    assert contract['nsvars']==27
    for key in ['H','alpha']:
        slots=contract['history_slots'][key]
        assert len(slots)==9 and len(set(slots))==9 and min(slots)>=1 and max(slots)<=27
    for name in ['design.md','research_log.md','STATUS.md','verification/RESULTS.md']:
        assert (ROOT/name).is_file()
    provenance=json.loads((ROOT/'verification'/'provenance.json').read_text())
    assert provenance['selected_source_sha256']==digest
    for name,old in provenance['host_sha256'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==old,name
    runs=[]
    for path in sorted((ROOT/'verification'/'runs').glob('selected_*')):
        status=json.loads((path/'run.json').read_text())
        manifest=json.loads((path/'input_manifest.json').read_text())
        assert manifest['source_sha256']==digest,path
        assert manifest['contract']['jprops']==[1]
        assert status['status']=='ok',path
        assert status['accepted_output_steps']==status['cycles']*status['steps_per_cycle']
        assert all(v is not None and v>=0. for v in status['history_min_change'].values())
        obs=np.genfromtxt(path/'observables.csv',delimiter=',',names=True)
        assert all(np.isfinite(obs[k]).all() for k in obs.dtype.names)
        assert (np.diff(obs['t'])>0).all()
        previous=None
        for f in sorted(path.glob('fields_cycle_*.npz')):
            a=np.load(f)
            assert np.isfinite(a['U']).all() and np.isfinite(a['svars']).all()
            d=a['U'][3::6]
            assert d.min()>-1e-8 and d.max()<1+1e-8,(f,d.min(),d.max())
            if previous is not None:
                for slots in contract['history_slots'].values():
                    idx=np.array(slots)-1
                    assert (a['svars'][:,idx]>=previous[:,idx]).all()
            previous=a['svars']
        runs.append({'name':path.name,'cycles':status['cycles'],'n':status['n'],
                     'steps_per_cycle':status['steps_per_cycle'],'status':'ok'})
    assert len(runs)==11,len(runs)
    for name in ['reject_before','reject_after','large_step_10','D0_statistics','constrained_notch','saturated_fatigue']:
        summary=json.loads((ROOT/'verification'/'global_checks'/name/'summary.json').read_text())
        assert summary['status']=='ok' and summary['source_sha256']==digest,name
        assert summary['history_min_change']['H']>=0 and summary['history_min_change']['alpha']>=0
    a=np.load(ROOT/'verification'/'global_checks'/'reject_before'/'final.npz')
    b=np.load(ROOT/'verification'/'global_checks'/'reject_after'/'final.npz')
    assert np.array_equal(a['U'],b['U']) and np.array_equal(a['svars'],b['svars'])
    for path in (ROOT/'verification'/'q2_reference').glob('n*_spc*/summary.json'):
        summary=json.loads(path.read_text());assert summary['status']=='ok',path
    result={'selected_source_sha256':digest,'host_hashes_unchanged':True,
            'selected_runs_finite_and_history_monotone':runs,
            'rollback_bitwise_equal':True,'required_global_checks_pass':True}
    write_json(ROOT/'verification'/'audit.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
