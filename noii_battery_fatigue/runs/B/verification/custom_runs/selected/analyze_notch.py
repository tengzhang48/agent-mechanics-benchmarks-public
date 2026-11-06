"""Read-only comparisons of frozen-source notched runs, with independent Q1 L2 norms."""
from pathlib import Path
import csv,json
import numpy as np
ROOT=Path(__file__).resolve().parent
CORNERS=np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])
def quadrature(X):
 for xi,eta in CORNERS/np.sqrt(3):
  s,t=CORNERS.T;N=(1+s*xi)*(1+t*eta)/4
  DN=np.array([s*(1+t*eta),t*(1+s*xi)]).T/4
  yield N,float(np.linalg.det(X.T@DN))
def main():
 rows=[];states=[];source=json.loads((ROOT/'source_bound_cases.json').read_text())['source_sha256'];names=['ux','uy','c','d','q','vfrac']
 for spc in [20,40,80,160,320,640]:
  folder=ROOT/('notch_distorted16_s3_c5_spc80' if spc==80 else f'notch_distorted16_s3_c2_spc{spc}')
  if not (folder/'run.json').exists():continue
  run=json.loads((folder/'run.json').read_text())
  if run['status']!='ok':rows.append(dict(spc=spc,status=run['status']));continue
  with np.load(folder/'fields_cycle_002.npz') as d:
   U=d['U'].reshape(-1,6).copy();sv=d['svars'].copy();nodes=d['nodes'].copy();elems=d['elems'].copy()
  obs=next(x for x in csv.DictReader((folder/'observables.csv').open()) if abs(float(x['t'])-16)<1e-10)
  rows.append(dict(spc=spc,status='ok',cutbacks=run['cutbacks'],alpha_max=float(sv[:,1::4].max()),H_max=float(sv[:,0::4].max()),intact_constraint_max=float(np.max(abs((U[:,3]-1)+U[:,5]))),**{k:float(v) for k,v in obs.items()}));states.append(dict(U=U,sv=sv,nodes=nodes,elems=elems))
 diffs=[]
 for i in range(len(states)-1):
  a,b=states[i:i+2];assert np.array_equal(a['nodes'],b['nodes']) and np.array_equal(a['elems'],b['elems']);du=a['U']-b['U'];l2=np.zeros(6);h2=np.zeros(4)
  for ei,en in enumerate(a['elems']):
   for gp,(N,w) in enumerate(quadrature(a['nodes'][en])):
    l2+=(N@du[en])**2*w;h2+=(a['sv'][ei].reshape(4,4)[gp]-b['sv'][ei].reshape(4,4)[gp])**2*w
  diffs.append(dict(coarse_spc=rows[i]['spc'],fine_spc=rows[i+1]['spc'],fields_max_abs={k:float(np.max(abs(du[:,j]))) for j,k in enumerate(names)},fields_l2={k:float(np.sqrt(l2[j])) for j,k in enumerate(names)},histories_gp_l2={k:float(np.sqrt(h2[j])) for j,k in enumerate(['H','alpha','W','psi'])},observable_abs={k:abs(rows[i][k]-rows[i+1][k]) for k in ['D_int','C_int','E_el','A_damage']}))
 orders=[]
 for i in range(len(diffs)-1):
  def order(a,b):return float(np.log2(a/b)) if a>0 and b>0 else None
  orders.append(dict(spc_triple=[rows[i]['spc'],rows[i+1]['spc'],rows[i+2]['spc']],fields_l2={k:order(diffs[i]['fields_l2'][k],diffs[i+1]['fields_l2'][k]) for k in names},field_d_max=order(diffs[i]['fields_max_abs']['d'],diffs[i+1]['fields_max_abs']['d']),observables={k:order(diffs[i]['observable_abs'][k],diffs[i+1]['observable_abs'][k]) for k in ['D_int','C_int','E_el','A_damage']}))
 result=dict(mesh='distorted',n=16,seed=3,cycles=2,source_sha256=source,runs=rows,successive_differences=diffs,observed_orders=orders,interpretation='Same-formulation temporal self-convergence; no independent notched PDE oracle is claimed.')
 (ROOT/'notch_time_refinement.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps(dict(runs=[dict(spc=r['spc'],D_int=r.get('D_int'),C_int=r.get('C_int'),E_el=r.get('E_el'),alpha_max=r.get('alpha_max')) for r in rows],differences=[dict(spc=[x['coarse_spc'],x['fine_spc']],d_L2=x['fields_l2']['d'],c_L2=x['fields_l2']['c'],d_max=x['fields_max_abs']['d'],D_int=x['observable_abs']['D_int']) for x in diffs]),indent=2))
if __name__=='__main__':main()
