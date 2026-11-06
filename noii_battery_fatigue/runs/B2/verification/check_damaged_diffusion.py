"""Resolve the low-mobility diffusion limit and diagnose consistent-mass ringing.

Uniform H=1e16, beta=0 and a top concentration step reduce the compiled Q4
problem to 1D diffusion with diffusivity kappa*D. The continuous-space BE
solution is a gamma-time average of the half-space erfc step solution.
The remote reflecting boundary has an exponentially negligible influence.
"""
from pathlib import Path
import sys
import json
import numpy as np
from scipy.linalg import solve_banded
from scipy.special import erfc,roots_genlaguerre,gammaln

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.common import compiled,properties
import fehost
import meshes


def discrete(n,diff,dt,steps,c0,amplitude):
    h=1/n
    md=np.full(n+1,2*h/3);md[[0,-1]]=h/3
    mo=h/6
    kd=np.full(n+1,2/h);kd[[0,-1]]=1/h
    ko=-1/h
    diag=md/dt+diff*kd;off=mo/dt+diff*ko
    band=np.zeros((3,n));band[1]=diag[:-1];band[0,1:]=off;band[2,:-1]=off
    u=np.full(n+1,c0)
    minimum=c0
    for _ in range(steps):
        rhs=md*u/dt
        rhs[1:]+=mo*u[:-1]/dt;rhs[:-1]+=mo*u[1:]/dt
        rhs[-2]-=off*(c0+amplitude)
        u[:-1]=solve_banded((1,1),band,rhs[:-1])
        u[-1]=c0+amplitude
        minimum=min(minimum,float(u.min()))
    return u,minimum


def host(n,c,dt,steps,H,sexact,amplitude):
    p=properties(c);nodes,elems=meshes.structured(n)
    el=fehost.ContractUEL(compiled(c),c,len(elems),dt);el.svars[:,0::3]=H
    el.svars[:,1::3]=.02  # exhausted-fatigue prehistory, currently unloaded
    v=np.zeros((len(nodes),6));v[:,2]=p['c0'];v[:,3]=1-sexact;v[:,5]=sexact
    top=np.where(nodes[:,1]==1)[0]
    n0=meshes.node_at(nodes,0,0);n1=meshes.node_at(nodes,1,0)
    bc={6*i+2:p['c0']+amplitude for i in top};bc.update({6*n0:0.,6*n0+1:0.,6*n1+1:0.})
    U,report=fehost.solve(nodes,elems,el,lambda t:bc,steps,dt,U0=v.ravel(),
        residual_scale=fehost.residual_weights(el,c))
    return U[2::6].reshape(n+1,n+1),report.summary()


def check_gamma_reference():
    # Dimensionless r=z/sqrt(D*dt). Check the BE resolvent recurrence
    # Y_N - d²Y_N/dr² = Y_(N-1) independently of the finite-element solve.
    r=np.r_[0.,np.geomspace(.01,150,80)]
    def values(steps,order):
        x,w=roots_genlaguerre(order,steps-1);w=w/w.sum()
        val=erfc(r[:,None]/(2*np.sqrt(x)))@w
        second=(r[:,None]/(2*np.sqrt(np.pi)*x**1.5)*np.exp(-r[:,None]**2/(4*x)))@w
        return val,second
    value,second=values(80,64);previous,_=values(79,64);finer,_=values(80,96)
    recurrence=float(np.max(abs(value-second-previous)))
    refinement=float(np.max(abs(value-finer)))
    assert recurrence<2e-12 and refinement<2e-12,(recurrence,refinement)
    return dict(BE_resolvent_recurrence_max_error=recurrence,quadrature_64_vs_96_max_error=refinement)


if __name__=='__main__':
    reference_checks=check_gamma_reference()
    c=fehost.load_contract(str(ROOT));p=properties(c);p['beta']=0.
    c={**c,'props':[p[k] for k in c['prop_names']]}
    H=1e16;sexact=1/(1+(1-p['kappa'])*H)
    diff=p['D']*(p['kappa']+(1-p['kappa'])*sexact**2)
    dt=.1;steps=80;end=dt*steps;amplitude=.45;c0=p['c0']
    gx,gw=roots_genlaguerre(64,steps-1);gw=gw/gw.sum()
    def be_exact(y):
        z=(1-np.asarray(y))[...,None]/(2*np.sqrt(diff*dt*gx))
        return c0+amplitude*np.sum(erfc(z)*gw,axis=-1)
    exact_mass=c0+amplitude*2*np.sqrt(diff*dt)/np.sqrt(np.pi)*np.exp(gammaln(steps+.5)-gammaln(steps))
    q,w=np.polynomial.legendre.leggauss(32)
    rows=[];host_rows=[]
    for n in (8,16,32,64,128,256,512,1024,2048):
        u,minimum=discrete(n,diff,dt,steps,c0,amplitude)
        y=(np.arange(n)[:,None]+(q[None]+1)/2)/n
        values=u[:-1,None]*(1-q[None])/2+u[1:,None]*(1+q[None])/2
        error=float(np.sqrt(np.sum((values-be_exact(y))**2*w[None]/(2*n))))
        mass=float((u.sum()-(u[0]+u[-1])/2)/n)
        rows.append(dict(n=n,c_min_over_steps=minimum,c_min_final=float(u.min()),
            C_int=mass,continuous_space_BE_C_int=float(exact_mass),L2_spatial_error=error))
        if n in (8,16):
            actual,report=host(n,c,dt,steps,H,sexact,amplitude)
            maximum=float(np.max(abs(actual-u[:,None])))
            assert maximum<2e-10,maximum
            host_rows.append(dict(n=n,compiled_vs_1D_max_error=maximum,solver=report))
    assert rows[0]['c_min_final']<c0-.01
    assert rows[-1]['L2_spatial_error']<rows[0]['L2_spatial_error']/100
    result=dict(status='diagnosed discretization limitation; compiled reduction verified',
        diffusivity=diff,time=end,diffusion_length=float(np.sqrt(diff*end)),
        boundary_layer_solution='BE resolvent = gamma-time mixture of erfc; 64 generalized Laguerre points',
        bottom_boundary_note='Distance one is over 298 diffusion lengths; reflecting-boundary correction is below roundoff',
        c0=c0,top_concentration=c0+amplitude,reference_checks=reference_checks,
        compiled_checks=host_rows,spatial_study=rows,
        finding='Consistent-mass Q1 diffusion can undershoot c0 when the low-mobility boundary layer is unresolved. Refinement resolves it. No positivity guarantee for concentration is claimed.')
    (ROOT/'verification/damaged_diffusion_resolution.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
