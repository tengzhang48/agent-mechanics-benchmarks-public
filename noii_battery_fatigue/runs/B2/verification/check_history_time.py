"""Analytical uniform-shear history/damage time-refinement oracle.

With beta=0, inactive fatigue, imposed gamma(t)=gamma_max*t and uniform
initial data, H(t)=H_final*t^2, d(t)=a*H/(1+a*H). The accumulated positive
energy has a separate analytical maximum. This checks the time-refined
continuum limit, including the declared one-increment history delay.
"""
from pathlib import Path
import sys
import json
import numpy as np
from scipy.optimize import brentq

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.common import compiled,properties
import fehost
import meshes


def run(steps):
    c=fehost.load_contract(str(ROOT));p=properties(c)
    gamma=.03;Hfinal=10.;a=1-p['kappa'];mu=p['E']/(2*(1+p['nu']))
    p.update(beta=0.,alpha_T=1.,G_c=p['ell']*mu*gamma**2/Hfinal)
    c={**c,'props':[p[k] for k in c['prop_names']]}
    nodes,elems=meshes.distorted(8,seed=6)
    dt=1/steps;el=fehost.ContractUEL(compiled(c),c,len(elems),dt)
    U0=np.zeros(len(nodes)*6);U0[2::6]=p['c0']
    def bc(t):
        vals={6*i:gamma*t*y for i,(x,y) in enumerate(nodes)}
        vals.update({6*i+1:0. for i in range(len(nodes))})
        return vals
    trajectory=[];minima=np.full(2,np.inf)
    def commit(X,conn,dpn,U,Uold):
        old=el.svars.copy();el.commit_step(X,conn,dpn,U,Uold)
        for i in range(2):
            change=el.svars[:,i::3]-old[:,i::3]
            minima[i]=min(minima[i],float(change.min()))
            assert np.all(change>=0)
    def output(step,t,U,metrics=None):
        trajectory.append([t,U[3::6].mean(),el.svars[:,0::3].mean(),el.svars[:,1::3].mean()])
    U,report=fehost.solve(nodes,elems,el,bc,steps,dt,U0=U0,
        residual_scale=fehost.residual_weights(el,c),commit_fn=commit,step_callback=output)
    exact_d=a*Hfinal/(1+a*Hfinal)
    split_H=Hfinal*(1-dt)**2
    split_d=a*split_H/(1+a*split_H)
    derror=float(np.max(abs(U[3::6]-exact_d)))
    lag_error=float(np.max(abs(U[3::6]-split_d)))
    assert lag_error<3e-10,lag_error
    herror=float(np.max(abs(el.svars[:,0::3]-Hfinal)))
    assert herror<2e-10,herror
    # W(H)=Gc/(2ell)*H*(kappa+a/(1+aH)^2), and it has a
    # single maximum on 0<=H<=10; alpha is that accumulated increase.
    derivative=lambda h:p['kappa']+a*(1-a*h)/(1+a*h)**3
    peak=brentq(derivative,0.,Hfinal)
    exact_alpha=p['G_c']/(2*p['ell'])*peak*(p['kappa']+a/(1+a*peak)**2)
    alpha_error=float(np.max(abs(el.svars[:,1::3]-exact_alpha)))
    out=ROOT/'verification/history_time';out.mkdir(exist_ok=True)
    np.savez_compressed(out/f'steps_{steps}.npz',trajectory=np.array(trajectory),
        columns=np.array(['t','d','H','alpha']),nodes=nodes,elems=elems,U=U,svars=el.svars)
    return dict(steps=steps,dt=dt,damage=float(U[3::6].mean()),continuum_damage=exact_d,
        damage_max_error=derror,analytical_split_damage_error=lag_error,H_max_error=herror,
        alpha=float(el.svars[:,1::3].mean()),continuum_alpha=exact_alpha,alpha_max_error=alpha_error,
        H_alpha_min_change=minima.tolist(),solver=report.summary())


if __name__=='__main__':
    rows=[run(n) for n in (16,32,64,128,256)]
    for old,new in zip(rows,rows[1:]):
        new['damage_rate']=float(np.log2(old['damage_max_error']/new['damage_max_error']))
        new['alpha_rate']=float(np.log2(old['alpha_max_error']/new['alpha_max_error']))
    assert rows[-1]['damage_max_error']<rows[0]['damage_max_error']/12
    assert rows[-1]['alpha_max_error']<rows[0]['alpha_max_error']/10
    result=dict(status='passed',scope='distorted n=8 seed=6, uniform prescribed shear, beta=0, inactive fatigue',runs=rows)
    (ROOT/'verification/history_time_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
