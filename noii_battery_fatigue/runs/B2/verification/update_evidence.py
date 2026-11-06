"""Refresh final-source audits and planned physical-space comparisons.

--watch waits for already-running studies. It launches no simulations and
never changes submitted sources, the host or the generator.
"""
from pathlib import Path
import sys
import json
import argparse
import time
import subprocess
from datetime import datetime,timezone

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.compare_fields import check,compare


def plans():
    rows=[]
    for n in (16,32,64):
        prod=f'verification/runs/D0_enhanced_accepted_H_structured_{n}_s0_spc80_cy20'
        for ref,tag in (('D0_Q4_structured_8_spc80_cy20','Q4n8'),
                        ('D0_Q3_structured_16_spc80_cy20','Q3n16')):
            rows.append((f'D0_n{n}_c20_{tag}',prod,'verification/reference_runs/'+ref,20))
        prod=f'verification/runs/D0_enhanced_accepted_H_distorted_{n}_s1_spc80_cy5'
        ref='verification/reference_runs/D0_Q4_distorted_8_spc80_cy5'
        rows.append((f'D0_dist{n}_c5_Q4n8',prod,ref,5))
    for n in (32,64):
        prod=f'verification/runs/V_severe_enhanced_accepted_H_structured_{n}_s0_spc160_cy2'
        ref='verification/reference_runs/V_severe_Q3_structured_16_spc160_cy2_acceptedH'
        rows.append((f'severe_n{n}_c2_dt160_Q3n16',prod,ref,2))
    rows.append(('D0_n32_c20_dt160_Q4n8',
        'verification/runs/D0_enhanced_accepted_H_structured_32_s0_spc160_cy20',
        'verification/reference_runs/D0_Q4_structured_8_spc160_cy20',20))
    rows.append(('severe_n32_c2_dt320_Q3n16',
        'verification/runs/V_severe_enhanced_accepted_H_structured_32_s0_spc320_cy2',
        'verification/reference_runs/V_severe_Q3_structured_16_spc320_cy2_acceptedH',2))
    rows.append(('severe_dist32_c2_dt160_Q3n16',
        'verification/runs/V_severe_enhanced_accepted_H_distorted_32_s3_spc160_cy2',
        'verification/reference_runs/V_severe_Q3_structured_16_spc160_cy2_acceptedH',2))
    rows.append(('severe_dist64_c2_dt160_Q3n16',
        'verification/runs/V_severe_enhanced_accepted_H_distorted_64_s3_spc160_cy2',
        'verification/reference_runs/V_severe_Q3_structured_16_spc160_cy2_acceptedH',2))
    for n in (16,32):
        rows.append((f'D0_n{n}_c20_dt40_Q4n8_acceptedH',
            f'verification/runs/D0_enhanced_accepted_H_structured_{n}_s0_spc40_cy20',
            'verification/reference_runs/D0_Q4_structured_8_spc40_cy20_acceptedH',20))
    for n in (32,64):
        for mesh,seed in (('structured',0),('distorted',3)):
            rows.append((f'severe_{mesh}{n}_c2_dt160_Q3n32',
                f'verification/runs/V_severe_enhanced_accepted_H_{mesh}_{n}_s{seed}_spc160_cy2',
                'verification/reference_runs/V_severe_Q3_structured_32_spc160_cy2_acceptedH',2))
            rows.append((f'severe_{mesh}{n}_c2_dt160_Q3n16q6',
                f'verification/runs/V_severe_enhanced_accepted_H_{mesh}_{n}_s{seed}_spc160_cy2',
                'verification/reference_runs/V_severe_Q3_structured_16_spc160_cy2_acceptedH_q6',2))
    for n in (16,32,64):
        for mesh,seed,cycles in (('structured',0,20),('distorted',1,5)):
            rows.append((f'D0_{mesh}{n}_c{cycles}_dt80_Q4n8_acceptedH',
                f'verification/runs/D0_enhanced_accepted_H_{mesh}_{n}_s{seed}_spc80_cy{cycles}',
                'verification/reference_runs/D0_Q4_structured_8_spc80_cy20_acceptedH',cycles))
    return rows


def main(watch=False,deadline='2026-10-21T05:00:00+00:00'):
    interpolation=check();out=ROOT/'verification/field_comparisons';out.mkdir(exist_ok=True)
    last=None
    while True:
        pending=[];changed=False
        for name,prod,ref,cycle in plans():
            target=out/(name+'.json')
            if target.exists() and json.loads(target.read_text())['comparison'].get('comparison_version')==2:continue
            required=[ROOT/prod/'run.json',ROOT/ref/'run.json']
            if not all(p.exists() and json.loads(p.read_text())['status']=='ok' for p in required):
                pending.append(name);continue
            result=dict(interpolation=interpolation,comparison=compare(ROOT/prod,ROOT/ref,cycle))
            target.write_text(json.dumps(result,indent=2)+'\n');changed=True
            print(datetime.now(timezone.utc).isoformat(),'Compared',name,flush=True)
        # Earlier manual comparisons retain useful Q3 mesh controls. Refresh
        # their observation diagnostics too, without deleting the old studies.
        for target in sorted(out.glob('*.json')):
            record=json.loads(target.read_text());rec=record['comparison']
            if rec.get('comparison_version')==2:continue
            result=dict(interpolation=interpolation,comparison=compare(
                ROOT/rec['production'],ROOT/rec['reference'],rec['cycle']))
            target.write_text(json.dumps(result,indent=2)+'\n');changed=True
            print(datetime.now(timezone.utc).isoformat(),'Updated',target.stem,flush=True)
        files=sorted((ROOT/'verification/runs').glob('*/run.json'))
        files+=sorted((ROOT/'verification/reference_runs').glob('*/run.json'))
        fingerprint=tuple((str(p),p.stat().st_mtime_ns) for p in files)
        if changed or fingerprint!=last:
            for script in ('audit_outputs.py','compare_enhancement.py','summarize.py'):
                log=ROOT/'tmp'/('update_'+script.replace('.py','.log'))
                with log.open('w') as stream:
                    subprocess.run([sys.executable,str(ROOT/'verification'/script)],cwd=ROOT,
                                   stdout=stream,stderr=subprocess.STDOUT,check=True)
            print(datetime.now(timezone.utc).isoformat(),'Refreshed evidence;',len(pending),'comparisons pending',flush=True)
            last=fingerprint
        if not watch or not pending or datetime.now(timezone.utc)>=datetime.fromisoformat(deadline):
            break
        time.sleep(20)
    print('Evidence updater finished; pending:',pending,flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--watch',action='store_true')
    ap.add_argument('--deadline',default='2026-10-21T05:00:00+00:00')
    main(**vars(ap.parse_args()))
