"""Higher-order structured 2D PDE comparison, using the unchanged generic host.

Reference-only nine-node Q2 macroelements. The scored element remains Q4.
"""
import argparse
import hashlib
import json
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from common import ROOT,fehost,load_element,raw_call,write_json


def variant():
    return dict(source='verification/q2_reference.f90',dof_per_node=4,
        dof_names=['ux','uy','c','d'],abaqus_dofs=[1,2,11,12],nsvars=75,
        history_slots={'H':list(range(1,76,3)),'alpha':list(range(2,76,3))},
        residual_scale={'ux':.1,'uy':.1,'c':.1,'d':1.},jprops=[])


def mesh(n):
    grid=np.linspace(0,1,2*n+1)
    xx,yy=np.meshgrid(grid,grid)
    nodes=np.column_stack([xx.ravel(),yy.ravel()])
    elems=np.array([[ (2*j+b)*(2*n+1)+2*i+a for b in range(3) for a in range(3)]
                    for j in range(n) for i in range(n)])
    return nodes,elems


def shape(xi,eta,h):
    lx=np.array([xi*(xi-1)/2,1-xi*xi,xi*(xi+1)/2])
    ly=np.array([eta*(eta-1)/2,1-eta*eta,eta*(eta+1)/2])
    dx=np.array([xi-.5,-2*xi,xi+.5]);dy=np.array([eta-.5,-2*eta,eta+.5])
    return np.outer(ly,lx).ravel(),np.column_stack([np.outer(ly,dx).ravel(),np.outer(dy,lx).ravel()])*2/h


def test():
    el,contract=load_element(1,.2,variant=variant());el.n_nodes=9;el.ndofel=36
    x,elems=mesh(1)
    rng=np.random.default_rng(18)
    uv=np.zeros((9,4));uv[:,:2]=rng.normal(0,.005,(9,2));uv[:,2]=.4+rng.random(9)*.1;uv[:,3]=.1+rng.random(9)*.1
    u=uv.ravel();old=u.copy();old[2::4]-=.03
    r,k,sv,pn=raw_call(el,x,u,old)
    fd=np.empty_like(k);eps=1e-7
    for j in range(36):
        du=np.zeros(36);du[j]=eps
        fd[:,j]=(raw_call(el,x,u+du,old)[0]-raw_call(el,x,u-du,old)[0])/(2*eps)
    error=float(np.linalg.norm(fd-k)/np.linalg.norm(k))
    assert error<1e-8,error
    patch=np.zeros((9,4));patch[:,2]=.35;patch[:,0]=.001*x[:,0]**2/2
    r,k,sv,pn=raw_call(el,x,patch.ravel(),patch.ravel())
    E,nu=10.,.28;mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu))
    flux=-.25*.35*2*(mu+lam)/3*.001
    expected=(np.outer([1/6,2/3,1/6],[-1.,0.,1.])*flux).ravel()
    patch_error=float(np.max(np.abs(r[2::4]-expected)))
    assert patch_error<1e-14,patch_error
    result={'Jacobian_relative_error':error,'quadratic_displacement_flux_residual_max_error':patch_error}
    write_json(ROOT/'verification'/'q2_reference'/'element_tests.json',result)
    print(json.dumps(result),flush=True)


def run(n,spc,cycles):
    prob=json.loads((ROOT/'host'/'dev_problem'/'d0.json').read_text());p=prob['parameters'];load=prob['loading']
    nodes,elems=mesh(n);dt=load['T']/spc
    el,contract=load_element(len(elems),dt,p,variant=variant());el.n_nodes=9;el.ndofel=36
    u0=np.zeros(4*len(nodes));u0[2::4]=p['c0']
    top=np.flatnonzero(np.isclose(nodes[:,1],1));origin=0;right=2*n
    def bc(frac):
        t=frac*cycles*load['T'];cb=p['c0']+load['A']*(1-np.cos(2*np.pi*t/load['T']))/2
        b={4*a+2:cb for a in top};b.update({4*origin:0.,4*origin+1:0.,4*right+1:0.});return b
    pts,ws=leggauss(5);geom=[(*shape(xi,eta,1/n),ws[i]*ws[j]/(4*n*n))
                    for j,eta in enumerate(pts) for i,xi in enumerate(pts)]
    mu=p['E']/(2*(1+p['nu']));lam=p['E']*p['nu']/((1+p['nu'])*(1-2*p['nu']))
    rows=[];history_min={'H':0.,'alpha':0.};previous=el.svars.copy()
    dest=ROOT/'verification'/'q2_reference'/f'n{n}_spc{spc}';dest.mkdir(parents=True,exist_ok=True)
    def output(step,frac,U,metrics=None):
        nonlocal previous
        uv=U.reshape(-1,4)[elems];total=np.zeros(3)
        for N,B,w in geom:
            val=np.einsum('a,maj->mj',N,uv);grad=np.einsum('ak,maj->mjk',B,uv[:,:,:2])
            ex=grad[:,0,0]-p['beta']*(val[:,2]-p['c0']);ey=grad[:,1,1]-p['beta']*(val[:,2]-p['c0'])
            es=(grad[:,0,1]+grad[:,1,0])/2
            g=(1-p['kappa'])*(1-val[:,3])**2+p['kappa']
            energy=g*(mu*(ex**2+ey**2+2*es**2)+lam*(ex+ey)**2/2)
            total+=w*np.array([val[:,3].sum(),val[:,2].sum(),energy.sum()])
        for name,sl in [('H',slice(0,None,3)),('alpha',slice(1,None,3))]:
            history_min[name]=min(history_min[name],float((el.svars[:,sl]-previous[:,sl]).min()))
        previous=el.svars.copy()
        rows.append([frac*cycles*load['T'],float(U[3::4].max()),*total])
        if step%spc==0:
            np.savez(dest/f'cycle{step//spc}.npz',nodes=nodes,elems=elems,U=U,svars=el.svars)
    start=time.time();status='ok';report=None
    try:
        U,report=fehost.solve(nodes,elems,el,bc,cycles*spc,dt,U0=u0,tol=1e-9,
            residual_scale=fehost.residual_weights(el,contract),step_callback=output)
    except Exception as exc:
        status=f'{type(exc).__name__}: {exc}'
    result={'n':n,'spc':spc,'cycles':cycles,'status':status,'wall_s':time.time()-start,
            'accepted_output_steps':len(rows),'history_min_change':history_min,
            'source_sha256':hashlib.sha256((ROOT/contract['source']).read_bytes()).hexdigest(),
            'final':rows[-1] if rows else None}
    if report is not None:result.update(report.summary())
    write_json(dest/'summary.json',result);np.savez(dest/'observables.npz',history=np.asarray(rows))
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=8);ap.add_argument('--spc',type=int,default=40)
    ap.add_argument('--cycles',type=int,default=5);ap.add_argument('--test',action='store_true');a=ap.parse_args()
    if a.test:test()
    else:run(a.n,a.spc,a.cycles)
