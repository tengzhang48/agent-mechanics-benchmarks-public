"""Analytical Q1 energy reconstruction on an exactly stress-free swelling field."""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'host'))
import meshes
from observables import observables


if __name__=='__main__':
    prm=json.loads((ROOT/'host/dev_problem/d0.json').read_text())['parameters']
    mu=prm['E']/(2*(1+prm['nu']))
    lam=prm['E']*prm['nu']/((1+prm['nu'])*(1-2*prm['nu']))
    rate=.1;scale=prm['beta']*rate;rows=[]
    for n in (8,16,32,64):
        nodes,elems=meshes.structured(n);x,y=nodes.T
        U=np.column_stack([scale*(x*x-y*y)/2,scale*x*y,prm['c0']+rate*x,np.zeros_like(x)])
        actual=observables(nodes,elems,U.ravel(),4,dict(ux=0,uy=1,c=2,d=3),prm)['E_el']
        analytical=scale**2*(3*mu+lam)/(24*n*n)
        assert abs(actual-analytical)<1e-12*analytical
        rows.append(dict(n=n,host_Q1_energy=float(actual),analytical_Q1_energy=analytical,
                         continuum_energy=0.,relative_formula_error=float(abs(actual-analytical)/analytical)))
    result=dict(status='passed',field='c=c0+a*x; ux=beta*a*(x*x-y*y)/2; uy=beta*a*x*y; d=0',
        explanation='The exact field is stress free. Sampling its exact nodal values and applying the host Q1 gradient creates O(h^2) energy.',runs=rows)
    (ROOT/'verification/energy_bridge_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
