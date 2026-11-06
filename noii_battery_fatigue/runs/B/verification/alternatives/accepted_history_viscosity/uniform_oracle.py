"""Every accepted compiled uniform state versus independent scalar bracketing."""
from pathlib import Path
import sys,json,numpy as np
from common import ROOT,compiled
from run_case import run
W=ROOT.parents[1];sys.path.insert(0,str(W/'verification'));from uniform_time_reference import solve
mc=compiled();rows=[]
for spc in [20,80,320]:
 out=ROOT/f'uniform_oracle_n8_c2_spc{spc}'
 r=run(W/'verification/problems/uniform_fatigue_floor.json',n=8,cycles=2,spc=spc,out=out,uniform_strain=.03,module_contract=mc)
 if r['status']!='ok':rows.append(dict(spc=spc,status=r['status']));continue
 data=np.genfromtxt(out/'increments.csv',delimiter=',',names=True)
 oracle=solve(spc,cycles=2,implicit='psilag',predictor_floor=.01,relaxation=.5,relaxation_cap=1000.,times=data['t'])
 ref=oracle['rows'];assert len(ref)==len(data)
 errors={}
 for name,col in [('H',1),('alpha',2),('W',3),('d',4)]:
  hi=data[name+'_max'];lo=data[name+'_min']
  errors[name]=dict(max_abs_error=float(max(np.max(abs(hi-ref[:,col])),np.max(abs(lo-ref[:,col])))),max_relative_error=float(max(np.max(abs(hi-ref[:,col])/np.maximum(abs(ref[:,col]),1e-20)),np.max(abs(lo-ref[:,col])/np.maximum(abs(ref[:,col]),1e-20)))),max_spatial_range=float(np.max(hi-lo)))
 rows.append(dict(spc=spc,status=r['status'],cutbacks=r['cutbacks'],accepted_subincrements=len(data),final_alpha=float(data['alpha_max'][-1]),scalar_final_alpha=oracle['final_alpha'],floor_onset_scalar=oracle['floor_onset'],errors=errors))
 (ROOT/'verification/uniform_scalar_comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(json.dumps(rows[-1]),flush=True)
