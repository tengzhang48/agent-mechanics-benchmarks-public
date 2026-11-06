"""Constitutive time benchmark with an independent scalar implicit reference.

Uniform prescribed c(t), ux=0, sigma_yy=0, homogeneous Neumann damage.
This isolates fatigue/history/damage time accuracy, without spatial error.
Concentration is prescribed everywhere; this is not a diffusion benchmark.
"""
import json
import time
import numpy as np
from scipy.optimize import brentq
from common import ROOT,fehost,load_element,write_json
from hostfe import CompiledElementState,solve_steps_staggered

P=dict(E=10.,nu=.28,beta=.035,c0=.35,kappa=1e-5,D=.14,chi=.25,
       ell=.12,G_c=3.5e-6,alpha_T=1e-6,k_fg=.5)
T=8.; AMP=.45; CYCLES=10
mu=P['E']/(2*(1+P['nu'])); lam=P['E']*P['nu']/((1+P['nu'])*(1-2*P['nu']))
ratio=lam/(2*mu+lam)
coef=.5*mu*P['beta']**2*(1+ratio)**2


def reference(spc):
    oldH=0.; olda=0.; oldW=0.; oldd=0.; result=[]
    for step in range(1,CYCLES*spc+1):
        t=step*T/spc; delta=AMP*(1-np.cos(2*np.pi*t/T))/2
        psi=coef*delta**2
        def state(d):
            g=(1-P['kappa'])*(1-d)**2+P['kappa']
            W=g*psi; alpha=olda+max(W-oldW,0.)
            factor=1. if alpha<=P['alpha_T'] else max(1e-15,max(0.,1-P['k_fg']*np.log10(alpha/P['alpha_T']))**2)
            H=max(oldH,2*P['ell']*psi/(P['G_c']*factor))
            A=(1-P['kappa'])*H
            return H,alpha,W,A/(1+A)
        def residual(d):
            return d-state(d)[3]
        # A [0,1] bracket also covers roundoff at the previous root on
        # unloading; using oldd as the lower endpoint can miss by a few ulps.
        d=brentq(residual,0.,1.,xtol=1e-14)
        oldH,olda,oldW,_=state(d); oldd=d
        result.append([t,d,oldH,olda,oldW])
    return np.asarray(result)


def uel(spc):
    x=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
    elems=np.array([[0,1,2,3]])
    el,contract=load_element(1,T/spc,P)
    u0=np.zeros(24); u0[2::6]=P['c0']; rows=[]
    def bc(frac):
        c=P['c0']+AMP*(1-np.cos(2*np.pi*frac*CYCLES))/2
        b={6*a:0. for a in range(4)}
        b.update({1:0.,7:0.})
        b.update({6*a+2:c for a in range(4)})
        return b
    def output(step,frac,u,metrics=None):
        rows.append([frac*CYCLES*T,float(u[3::6].mean()),
                     float(el.svars[0,0::3].mean()),float(el.svars[0,1::3].mean()),
                     float(el.svars[0,2::3].mean())])
    u,report=solve_steps_staggered(x,elems,6,el.element_rk,bc,[list(range(6))],CYCLES*spc,
        dt=T/spc,state=CompiledElementState([el]),batch_fn=el.element_rk_batch,
        commit_fn=el.commit_step,U0=u0,tol=1e-11,max_field_newton=25,max_cutbacks=10,
        residual_scale=fehost.residual_weights(el,contract),step_callback=output,return_report=True)
    return np.asarray(rows),report.summary()


def main():
    start=time.time(); ref=reference(5120)
    dest=ROOT/'verification'/'uniform_history'; dest.mkdir(parents=True,exist_ok=True)
    np.savez(dest/'reference.npz',history=ref,parameters=np.array(json.dumps(P)))
    result=[]
    for spc in [40,80,160,320,640]:
        values,stats=uel(spc)
        exact=np.column_stack([np.interp(values[:,0],ref[:,0],ref[:,j]) for j in range(1,5)])
        errors=np.max(np.abs(values[:,1:]-exact),axis=0)
        row={'spc':spc,'cycles':CYCLES,'final':values[-1].tolist(),
             'max_abs_error':dict(zip(['d','H','alpha','W'],errors.tolist())),**stats}
        result.append(row)
        np.savez(dest/f'uel_spc{spc}.npz',history=values)
        print(json.dumps(row),flush=True)
    # Independently check the reference temporal resolution.
    ref2=reference(2560)
    comparison=ref[1::2,1:]-ref2[:,1:]
    summary={'parameters':P,'reference_spc':5120,'cycles':CYCLES,
             'reference_2560_5120_max_change':dict(zip(['d','H','alpha','W'],np.max(np.abs(comparison),axis=0).tolist())),
             'reference_final':ref[-1].tolist(),'runs':result,'wall_s':time.time()-start}
    write_json(dest/'summary.json',summary)


if __name__=='__main__':
    main()
