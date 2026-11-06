"""Independent smooth 2D manufactured PDE, with analytical body sources.

Sympy differentiates the stated continuum laws, not the element residual.
The unchanged host solves the resulting monolithic problem. A frozen spatial
history isolates spatial consistency. Both auxiliary fields have no BCs.
"""
import argparse
import json
import time
import hashlib

import numpy as np
import sympy as sy
from numpy.polynomial.legendre import leggauss

from common import ROOT, fehost, meshes, load_element, write_json
from hostfe import CompiledElementState, solve_steps_staggered


def exact_functions():
    x, y = sy.symbols('x y', real=True)
    pi = sy.pi
    p = dict(E=10., nu=.28, beta=.035, c0=.35, kappa=1e-5, D=.14,
             chi=.25, ell=.12, G_c=1., alpha_T=1e4, k_fg=.5)
    mu = p['E']/(2*(1+p['nu']))
    lam = p['E']*p['nu']/((1+p['nu'])*(1-2*p['nu']))
    ux = .005*sy.sin(pi*x)*sy.sin(pi*y)+.002*x
    uy = .004*sy.cos(pi*x)*sy.sin(pi*y)-.001*y
    c = .43+.06*sy.sin(pi*x)*sy.sin(pi*y)+.03*x
    d = .25+.1*sy.sin(pi*x)*sy.sin(pi*y)
    h = .7+.2*sy.sin(2*pi*x)*sy.cos(pi*y)
    g = (1-p['kappa'])*(1-d)**2+p['kappa']
    ex = sy.diff(ux,x)-p['beta']*(c-p['c0'])
    ey = sy.diff(uy,y)-p['beta']*(c-p['c0'])
    es = (sy.diff(ux,y)+sy.diff(uy,x))/2
    sx, syy, sxy = g*(2*mu*ex+lam*(ex+ey)), g*(2*mu*ey+lam*(ex+ey)), g*2*mu*es
    q = g*(2*(mu+lam)/3)*(ex+ey)
    z = (1-p['kappa'])*h*(1-d)
    jx, jy = g*p['D']*sy.diff(c,x)-p['chi']*c*sy.diff(q,x), g*p['D']*sy.diff(c,y)-p['chi']*c*sy.diff(q,y)
    f = [-sy.diff(sx,x)-sy.diff(sxy,y), -sy.diff(sxy,x)-sy.diff(syy,y),
         -sy.diff(jx,x)-sy.diff(jy,y),
         -p['ell']**2*(sy.diff(d,x,2)+sy.diff(d,y,2))+d-z]
    expressions = [ux,uy,c,d,q,z,h,*f,sy.diff(q,x),sy.diff(q,y),jx,jy]
    fn = sy.lambdify((x,y),expressions,modules='numpy',cse=True)

    def evaluate(points):
        shape=points.shape[:-1]
        return np.stack([np.broadcast_to(a,shape) for a in
                         fn(points[...,0],points[...,1])],axis=-1)
    return p,evaluate


def geometry(nodes,elems,order):
    X=nodes[elems]
    pts,weights=leggauss(order)
    s=np.array([-1,1,1,-1]); t=np.array([-1,-1,1,1])
    for iy,eta in enumerate(pts):
        for ix,xi in enumerate(pts):
            N=.25*(1+s*xi)*(1+t*eta)
            dn=np.column_stack([.25*s*(1+t*eta),.25*t*(1+s*xi)])
            J=np.einsum('mai,aj->mij',X,dn)
            B=np.einsum('ai,mij->maj',dn,np.linalg.inv(J))
            yield N,B,np.linalg.det(J)*weights[ix]*weights[iy],np.einsum('a,mai->mi',N,X)


def run(n, mesh, lump, pressure, evaluate, p):
    nodes,elems=meshes.structured(n) if mesh=='structured' else meshes.distorted(n,seed=2)
    el,contract=load_element(len(elems),dt=1.,parameters=p)
    el._jprops=np.array([lump,pressure],dtype=np.int32)
    u0=evaluate(nodes)[:,:6].ravel()
    for ip,(_,_,_,points) in enumerate(geometry(nodes,elems,3)):
        el.svars[:,3*ip]=evaluate(points)[:,6]
        el.svars[:,3*ip+2]=1.
    boundary=meshes.boundary_nodes(nodes)
    bc={6*a+j:u0[6*a+j] for a in boundary for j in range(4)}
    force=np.zeros((len(elems),4,6))
    for N,B,w,points in geometry(nodes,elems,5):
        f=evaluate(points)[:,7:11]
        force[:,:,:4]+=N[None,:,None]*w[:,None,None]*f[:,None,:]
    force=force.reshape(len(elems),24)
    def loaded_batch(coords,u,du):
        r,k=el.element_rk_batch(coords,u,du)
        return r-force,k
    start=time.time()
    U,report=solve_steps_staggered(nodes,elems,6,el.element_rk,lambda frac:bc,
        [list(range(6))],1,batch_fn=loaded_batch,commit_fn=el.commit_step,
        state=CompiledElementState([el]),dt=1.,max_cutbacks=10,
        residual_scale=fehost.residual_weights(el,contract),U0=u0,tol=1e-10,
        max_field_newton=25,return_report=True)
    uv=U.reshape(-1,6)[elems]
    error=np.zeros(6); gradqerr=0.; fluxerr=0.
    for N,B,w,points in geometry(nodes,elems,5):
        exact=evaluate(points)
        val=np.einsum('a,maj->mj',N,uv)
        gradq=np.einsum('mak,ma->mk',B,uv[:,:,4])
        gradc=np.einsum('mak,ma->mk',B,uv[:,:,2])
        g=(1-p['kappa'])*(1-val[:,3])**2+p['kappa']
        flux=p['D']*g[:,None]*gradc-p['chi']*val[:,2,None]*gradq
        error+=np.einsum('m,mj->j',w,(val-exact[:,:6])**2)
        gradqerr+=float(np.sum(w[:,None]*(gradq-exact[:,11:13])**2))
        fluxerr+=float(np.sum(w[:,None]*(flux-exact[:,13:15])**2))
    result={'n':n,'mesh':mesh,'lumping':lump,'mixed_pressure':pressure,
            'source_sha256':hashlib.sha256((ROOT/'element'/'uel.f90').read_bytes()).hexdigest(),
            'L2':dict(zip(contract['dof_names'],np.sqrt(error).tolist())),
            'grad_q_L2':gradqerr**.5,'flux_L2':fluxerr**.5,
            'wall_s':round(time.time()-start,3),**report.summary()}
    path=ROOT/'verification'/'manufactured'/f'{mesh}_n{n}_lump{lump}_p{pressure}.json'
    write_json(path,result)
    print(json.dumps(result),flush=True)
    return result


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--sizes',type=int,nargs='+',default=[8,16,32])
    ap.add_argument('--meshes',nargs='+',default=['structured','distorted'])
    ap.add_argument('--lumping',type=int,default=1)
    ap.add_argument('--pressure',type=int,default=0)
    a=ap.parse_args()
    p,ev=exact_functions()
    for mesh in a.meshes:
        for n in a.sizes:
            run(n,mesh,a.lumping,a.pressure,ev,p)


if __name__=='__main__':
    main()
