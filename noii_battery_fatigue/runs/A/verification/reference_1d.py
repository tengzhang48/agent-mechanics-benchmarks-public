"""Independent cell-centred finite-volume reference and a host strip test.

ux=0 throughout a narrow strip, uy=0 on its base. All fields are uniform in x.
Equilibrium gives sigma_yy=0 and eliminates uy exactly in the reference.
This tests stress-assisted diffusion, swelling, fatigue and damage together.
The FV reference uses no element source, auxiliary q or auxiliary z.
"""
import argparse
import json
import time
import hashlib

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

from common import ROOT, fehost, load_element, write_json
from hostfe import CompiledElementState, solve_steps_staggered
from observables import observables

P=dict(E=10.,nu=.28,beta=.035,c0=.35,kappa=1e-5,D=.14,chi=.25,
       ell=.12,G_c=1e-4,alpha_T=5e-4,k_fg=.5)
PERIOD=8.; AMP=.3


def history(c,d,sv):
    p=P; mu=p['E']/(2*(1+p['nu'])); lam=p['E']*p['nu']/((1+p['nu'])*(1-2*p['nu']))
    ratio=lam/(2*mu+lam); delta=c-p['c0']; bulk=lam+2*mu/3
    tr=p['beta']*(ratio-1)*delta; dev=-p['beta']*(1+ratio)*delta
    trc=p['beta']*(ratio-1); devc=-p['beta']*(1+ratio)
    psi=.5*bulk*np.maximum(tr,0)**2+.5*mu*dev**2
    pc=bulk*np.maximum(tr,0)*trc+mu*dev*devc
    g=(1-p['kappa'])*(1-d)**2+p['kappa']; gd=-2*(1-p['kappa'])*(1-d)
    w=g*psi
    alpha=sv[:,1]+np.maximum(w-sv[:,2],0)
    ac=np.where(w>sv[:,2],g*pc,0.); ad=np.where(w>sv[:,2],gd*psi,0.)
    bracket=np.ones_like(c); mask=alpha>p['alpha_T']
    bracket[mask]=np.maximum(0,1-p['k_fg']*np.log10(alpha[mask]/p['alpha_T']))
    f=np.maximum(1e-15,bracket**2)
    fa=np.zeros_like(c); active=mask&(bracket**2>1e-15)
    fa[active]=-2*bracket[active]*p['k_fg']/(alpha[active]*np.log(10.))
    trial=2*p['ell']*psi/(p['G_c']*f)
    hc=2*p['ell']*pc/(p['G_c']*f)-trial*fa/f*ac
    hd=-trial*fa/f*ad
    hc=np.where(trial>sv[:,0],hc,0.); hd=np.where(trial>sv[:,0],hd,0.)
    h=np.maximum(trial,sv[:,0]); A=(1-p['kappa'])*h; ih=1/(1+A)
    return np.column_stack([h,alpha,w]),g,gd,A,ih,(1-p['kappa'])*hc,(1-p['kappa'])*hd


def fv_evaluate(v,old,sv,dt,cb,L,lag=False):
    n=len(v)//2; h=1/n; c=v[:n]; d=v[n:]; p=P
    sn,g,gd,A,ih,Ac,Ad=history(c,d,sv)
    if lag:
        A=(1-p['kappa'])*sv[:,0]
        ih=1/(1+A)
        Ac=np.zeros(n); Ad=np.zeros(n)
    mu=p['E']/(2*(1+p['nu'])); lam=p['E']*p['nu']/((1+p['nu'])*(1-2*p['nu']))
    pc=4*mu*(mu+lam)*p['beta']/(3*(2*mu+lam))
    q=-pc*g*(c-p['c0']); qc=-pc*g; qd=-pc*gd*(c-p['c0'])
    dc=(c[1:]-c[:-1])/h; dq=(q[1:]-q[:-1])/h
    ga=.5*(g[1:]+g[:-1]); ca=.5*(c[1:]+c[:-1])
    flux=p['D']*ga*dc-p['chi']*ca*dq
    flc=-p['D']*ga/h-p['chi']*.5*dq+p['chi']*ca*qc[:-1]/h
    frc=p['D']*ga/h-p['chi']*.5*dq-p['chi']*ca*qc[1:]/h
    fld=.5*p['D']*gd[:-1]*dc+p['chi']*ca*qd[:-1]/h
    frd=.5*p['D']*gd[1:]*dc-p['chi']*ca*qd[1:]/h
    top=g[-1]*(p['D']+p['chi']*pc*.5*(c[-1]+cb))*(cb-c[-1])/(h/2)
    tc=g[-1]*(p['chi']*pc*.5*(cb-c[-1])-(p['D']+p['chi']*pc*.5*(c[-1]+cb)))/(h/2)
    td=gd[-1]*(p['D']+p['chi']*pc*.5*(c[-1]+cb))*(cb-c[-1])/(h/2)
    fluxes=np.r_[0.,flux,top]
    rc=(c-old[:n])/dt-np.diff(fluxes)/h
    lap=L@d
    rd=d-A*ih-p['ell']**2*ih*lap
    # Each interface contributes opposite signed flux to adjacent volumes.
    row=np.repeat(np.arange(n-1),2); row=np.column_stack([np.arange(n-1),np.arange(1,n)]).ravel()
    colleft=np.repeat(np.arange(n-1),2); colright=colleft+1
    sign=np.tile([-1.,1.],n-1)/h
    jcc=sp.coo_matrix((np.r_[sign*np.repeat(flc,2),sign*np.repeat(frc,2),[-tc/h]],
                      (np.r_[row,row,n-1],np.r_[colleft,colright,n-1])),shape=(n,n)).tocsr()+sp.eye(n)/dt
    jcd=sp.coo_matrix((np.r_[sign*np.repeat(fld,2),sign*np.repeat(frd,2),[-td/h]],
                      (np.r_[row,row,n-1],np.r_[colleft,colright,n-1])),shape=(n,n)).tocsr()
    coeff=(p['ell']**2*lap-1)*ih**2
    jdc=sp.diags(coeff*Ac)
    jdd=sp.eye(n)-p['ell']**2*sp.diags(ih)@L+sp.diags(coeff*Ad)
    return np.r_[rc,rd],sp.bmat([[jcc,jcd],[jdc,jdd]],format='csr'),sn


def fv(n,spc,cycles,lag=False):
    start=time.time(); dt=PERIOD/spc; h=1/n
    diag=np.full(n,-2.); diag[[0,-1]]=-1.
    L=sp.diags([np.ones(n-1),diag,np.ones(n-1)],[-1,0,1],format='csr')/h**2
    v=np.r_[np.full(n,P['c0']),np.zeros(n)]; sv=np.zeros((n,3)); rows=[]; cutbacks=0; maxit=0
    previous_t=0.
    def increment(v0,s0,t,stepdt):
        cb=P['c0']+AMP*(1-np.cos(2*np.pi*t/PERIOD))/2
        vv=v0.copy()
        for it in range(35):
            r,k,sn=fv_evaluate(vv,v0,s0,stepdt,cb,L,lag)
            norm=np.linalg.norm(r,np.inf)
            if norm<2e-10:
                return vv,sn,it
            update=spsolve(k,-r); factor=1.
            for ls in range(25):
                trial=vv+factor*update
                rt,_,_=fv_evaluate(trial,v0,s0,stepdt,cb,L,lag)
                if np.linalg.norm(rt,np.inf)<norm:
                    vv=trial
                    break
                factor*=.5
            else:
                raise RuntimeError('FV line search')
        raise RuntimeError('FV Newton')
    def advance(v,s,t0,t1,depth=0):
        nonlocal cutbacks,maxit
        try:
            vv,ss,it=increment(v,s,t1,t1-t0)
            maxit=max(maxit,it)
            assert (ss[:,:2]>=s[:,:2]).all()
            return vv,ss
        except RuntimeError:
            if depth>=10:
                raise
            cutbacks+=1
            mid=.5*(t0+t1)
            vv,ss=advance(v,s,t0,mid,depth+1)
            return advance(vv,ss,mid,t1,depth+1)
    mu=P['E']/(2*(1+P['nu'])); lam=P['E']*P['nu']/((1+P['nu'])*(1-2*P['nu']))
    ratio=lam/(2*mu+lam)
    for step in range(1,cycles*spc+1):
        t=step*dt; v,sv=advance(v,sv,previous_t,t); previous_t=t
        c=v[:n]; d=v[n:]; g=(1-P['kappa'])*(1-d)**2+P['kappa']
        ex=-P['beta']*(c-P['c0']); ey=-ratio*ex
        e=g*(mu*(ex**2+ey**2)+.5*lam*(ex+ey)**2)
        rows.append([t,float(d.max()),float(d.mean()),float(c.mean()),float(e.mean())])
    result={'method':'lag_independent_FV' if lag else 'independent_FV','n':n,'spc':spc,'cycles':cycles,'wall_s':time.time()-start,
            'cutbacks':cutbacks,'max_newton_updates':maxit,'final':rows[-1],
            'parameters':P,'period':PERIOD,'amplitude':AMP}
    prefix='lag_fv' if lag else 'fv'
    dest=ROOT/'verification'/'reference_1d'/f'{prefix}_n{n}_spc{spc}'
    write_json(dest/'summary.json',result)
    np.savez(dest/'result.npz',history=np.asarray(rows),y=(np.arange(n)+.5)/n,c=v[:n],d=v[n:],svars=sv)
    print(json.dumps(result),flush=True)


def strip(n,spc,cycles,lag=False):
    start=time.time(); width=1/n
    y=np.linspace(0,1,n+1)
    nodes=np.column_stack([np.tile([0.,width],n+1),np.repeat(y,2)])
    elems=np.array([[2*j,2*j+1,2*j+3,2*j+2] for j in range(n)])
    variant={'source':'verification/alternatives/uel_lagged.f90'} if lag else None
    dt=PERIOD/spc; el,contract=load_element(n,dt,P,variant=variant)
    u0=np.zeros(6*len(nodes)); u0[2::6]=P['c0']
    def bc(frac):
        t=frac*cycles*PERIOD; cb=P['c0']+AMP*(1-np.cos(2*np.pi*t/PERIOD))/2
        result={6*a:0. for a in range(len(nodes))}
        result.update({1:0.,7:0.,6*(len(nodes)-2)+2:cb,6*(len(nodes)-1)+2:cb})
        return result
    rows=[]
    def output(step,frac,U,metrics=None):
        ob=observables(nodes,elems,U,6,{'ux':0,'uy':1,'c':2,'d':3},P)
        rows.append([frac*cycles*PERIOD,float(U[3::6].max()),
                     ob['D_int']/width,ob['C_int']/width,ob['E_el']/width])
    U,report=solve_steps_staggered(nodes,elems,6,el.element_rk,bc,[list(range(6))],cycles*spc,
          dt=dt,state=CompiledElementState([el]),batch_fn=el.element_rk_batch,
          commit_fn=el.commit_step,U0=u0,tol=1e-10,max_field_newton=25,max_cutbacks=10,
          residual_scale=fehost.residual_weights(el,contract),step_callback=output,return_report=True)
    result={'method':'lag_UEL_strip' if lag else 'UEL_strip','n':n,'spc':spc,'cycles':cycles,'wall_s':time.time()-start,
            'source_sha256':hashlib.sha256((ROOT/contract['source']).read_bytes()).hexdigest(),
            'final':rows[-1],**report.summary()}
    prefix='lag_uel' if lag else 'uel'
    dest=ROOT/'verification'/'reference_1d'/f'{prefix}_n{n}_spc{spc}'
    write_json(dest/'summary.json',result)
    uv=U.reshape(-1,6).reshape(n+1,2,6).mean(axis=1)
    np.savez(dest/'result.npz',history=np.asarray(rows),y=y,c=uv[:,2],d=uv[:,3],svars=el.svars,U=U)
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('method',choices=['fv','uel'])
    ap.add_argument('--n',type=int,default=128)
    ap.add_argument('--spc',type=int,default=80)
    ap.add_argument('--cycles',type=int,default=5)
    ap.add_argument('--lag',action='store_true')
    a=ap.parse_args()
    if a.method=='fv':
        fv(a.n,a.spc,a.cycles,a.lag)
    else:
        strip(a.n,a.spc,a.cycles,a.lag)
