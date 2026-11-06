"""Print only new milestones from already-running studies; no simulation work."""
from pathlib import Path
import json
import subprocess

ROOT=Path(__file__).resolve().parents[1]
NAMES=('final_enhanced64_long','final_enhanced32_long160','D0_enhanced_long320',
       'dist64_enhanced_final','severe64_enhanced_final160','severe_dist64_enhanced_final160',
       'reference_Q3n16_long80','reference_Q4_long160','reference_severe16_320',
       'reference_Q4_long40_acceptedH','reference_Q4_long80_acceptedH',
       'reference_severe32_160','reference_severe16_160_q6','evidence_updater')


if __name__=='__main__':
    cache=ROOT/'tmp/progress_seen.json'
    seen=json.loads(cache.read_text()) if cache.exists() else {}
    print(subprocess.check_output(['date','-u','+%Y-%m-%dT%H:%M:%SZ'],text=True).strip())
    changed=False
    for name in NAMES:
        path=ROOT/'tmp'/(name+'.log')
        lines=path.read_text().splitlines() if path.exists() else []
        current=lines[-1] if lines else 'No cycle saved yet'
        if seen.get(name)==current:continue
        seen[name]=current;changed=True
        if current.startswith('{'):
            try:
                r=json.loads(current)
                current=json.dumps({k:r[k] for k in ('status','n','cycles','steps_per_cycle','spc','accepted_output_steps','cutbacks','d_max_over_run') if k in r})
            except json.JSONDecodeError:pass
        print(name+': '+current)
    result=ROOT/'verification/results.json'
    if result.exists():
        r=json.loads(result.read_text())
        counts=[len(r['completed_final_source_runs']),len(r['pending_final_source_runs']),len(r['failed_final_source_runs'])]
        if seen.get('counts')!=counts:
            print('Final-source completed/pending/failed:',counts);seen['counts']=counts;changed=True
    if not changed:print('No new milestone; studies continue.')
    cache.write_text(json.dumps(seen,indent=2)+'\n')
