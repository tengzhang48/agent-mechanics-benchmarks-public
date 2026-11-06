"""Compare saved Q1 primary fields with the independent Qp reference.

Reference values are interpolated in physical space using an inverse bilinear
geometry map, including on distorted meshes. No UEL or model declaration is
used. A polynomial reproduction check validates the interpolation bridge.
"""
from pathlib import Path
import sys
import json
import argparse
import hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.reference import shape,enhanced_strain,GP
from verification.high_order_reference import mesh_qp


def interpolate(nodes,elems,values,points,n,order):
    """Vectorized point location followed by tensor Lagrange interpolation."""
    nn=(order+1)**2
    corners=elems[:,[0,order,nn-1,order*(order+1)]]
    X=nodes[corners]
    base=np.minimum(np.maximum(np.floor(points*n).astype(int),0),n-1)
    chosen=np.full(len(points),-1,int);natural=np.zeros((len(points),2))
    signs=np.array([[-1.,-1.],[1.,-1.],[1.,1.],[-1.,1.]])
    for dx,dy in [(0,0)]+[(i,j) for j in (-1,0,1) for i in (-1,0,1) if i or j]:
        idx=np.where(chosen<0)[0]
        if not len(idx):break
        grid=base[idx]+[dx,dy]
        valid=np.all((grid>=0)&(grid<n),axis=1)
        idx=idx[valid];grid=grid[valid]
        eid=grid[:,1]*n+grid[:,0]
        geom=X[eid]
        rs=np.zeros((len(idx),2))
        for _ in range(12):
            N=(1+rs[:,0,None]*signs[None,:,0])*(1+rs[:,1,None]*signs[None,:,1])/4
            dn=np.empty((len(idx),4,2))
            dn[:,:,0]=signs[None,:,0]*(1+rs[:,1,None]*signs[None,:,1])/4
            dn[:,:,1]=signs[None,:,1]*(1+rs[:,0,None]*signs[None,:,0])/4
            residual=np.einsum('pa,pai->pi',N,geom)-points[idx]
            J=np.einsum('pai,par->pir',geom,dn)
            det=J[:,0,0]*J[:,1,1]-J[:,0,1]*J[:,1,0]
            safe=np.where(abs(det)>1e-15,det,1e-15)
            step=np.column_stack((J[:,1,1]*residual[:,0]-J[:,0,1]*residual[:,1],
                                  -J[:,1,0]*residual[:,0]+J[:,0,0]*residual[:,1]))/safe[:,None]
            rs=np.clip(rs-np.clip(step,-2,2),-3,3)
            if np.max(abs(step),initial=0)<2e-13:break
        # Wrong candidate cells need not admit a well-behaved inverse outside
        # their domain. Accept only a converged inverse lying inside the cell.
        inside=(np.max(abs(rs),axis=1)<=1+2e-10)&(np.linalg.norm(residual,axis=1)<1e-10)
        chosen[idx[inside]]=eid[inside];natural[idx[inside]]=rs[inside]
    assert np.all(chosen>=0),points[chosen<0]
    pts=np.linspace(-1,1,order+1)
    lag=np.ones((len(points),2,order+1))
    for i in range(order+1):
        for j in range(order+1):
            if i!=j:lag[:,:,i]*=(natural-pts[j])/(pts[i]-pts[j])
    N=(lag[:,1,:,None]*lag[:,0,None,:]).reshape(len(points),nn)
    return np.einsum('pa,paj->pj',N,values[elems[chosen]])


def check():
    errors={}
    for distorted in (False,True):
        for order in (2,3,4):
            nodes,elems=mesh_qp(8,order,distorted,seed=2)
            def exact(x):
                a,b=x.T
                return np.column_stack([1+a*a+.3*a*b, a-b*b+.2*a*b, .35+.02*a])
            points=np.random.default_rng(731).uniform(0,1,(250,2))
            val=interpolate(nodes,elems,exact(nodes),points,8,order)
            err=float(np.max(abs(val-exact(points))))
            assert err<3e-13,err
            errors[f"Q{order}_{'distorted' if distorted else 'structured'}"]=err
    return dict(status="passed",polynomial_max_abs_errors=errors)


def compare(production,reference,cycle):
    p=Path(production);r=Path(reference)
    ps=json.loads((p/'run.json').read_text());rs=json.loads((r/'run.json').read_text())
    assert ps['status']==rs['status']=='ok'
    assert cycle<=min(ps['cycles'],rs['cycles'])
    a=np.load(p/f'fields_cycle_{cycle:03d}.npz',allow_pickle=False)
    b=np.load(r/f'fields_cycle_{cycle:03d}.npz',allow_pickle=False)
    names=a['dof_names'].tolist();idx=[names.index(k) for k in ('ux','uy','c','d')]
    av=a['U'].reshape(-1,len(names));bv=b['U'].reshape(-1,4)
    X=a['nodes'][a['elems']];vals=av[a['elems']]
    error=np.zeros(4);norm=np.zeros(4)
    q,w=np.polynomial.legendre.leggauss(4)
    for i,xi in enumerate(q):
        for j,eta in enumerate(q):
            N,dn=shape(xi,eta)
            weights=np.linalg.det(np.einsum('eai,ar->eir',X,dn))*w[i]*w[j]
            points=np.einsum('eai,a->ei',X,N)
            actual=np.einsum('eaj,a->ej',vals[:,:,idx],N)
            target=interpolate(b['nodes'],b['elems'],bv,points,rs['n'],rs['order'])
            error+=np.einsum('e,ej->j',weights,(actual-target)**2)
            norm+=np.einsum('e,ej->j',weights,target**2)
    # Independent energy recovery, complementary to the host's Q1 observable.
    prm=ps['parameters'];E=prm['E'];nu=prm['nu']
    mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu))
    sampled=interpolate(b['nodes'],b['elems'],bv,a['nodes'],rs['n'],rs['order'])
    reference_Q1_energy=0.
    for xi,eta in GP:
        N,dn=shape(xi,eta)
        jac=np.einsum('eai,ar->eir',X,dn);B=np.einsum('ar,eri->eai',dn,np.linalg.inv(jac))
        v=sampled[a['elems']]
        gu=np.einsum('eai,eaj->eij',v[:,:,:2],B)
        cg=v[:,:,2]@N;dg=v[:,:,3]@N
        ee=(gu+np.swapaxes(gu,1,2))/2-prm['beta']*(cg-prm['c0'])[:,None,None]*np.eye(2)
        g=(1-prm['kappa'])*(1-dg)**2+prm['kappa']
        psi=mu*np.sum(ee*ee,axis=(1,2))+lam*np.trace(ee,axis1=1,axis2=2)**2/2
        reference_Q1_energy+=float(np.sum(np.linalg.det(jac)*g*psi))
    enhanced_energy=0.
    for coords,v in zip(X,vals):
        correction=enhanced_strain(coords,v,prm)
        for k,(xi,eta) in enumerate(GP):
            N,dn=shape(xi,eta);J=coords.T@dn;B=dn@np.linalg.inv(J)
            gu=v[:,:2].T@B+correction[k]
            ee=(gu+gu.T)/2-prm['beta']*(N@v[:,2]-prm['c0'])*np.eye(2)
            g=(1-prm['kappa'])*(N@v[:,5])**2+prm['kappa']
            enhanced_energy+=np.linalg.det(J)*g*(mu*np.sum(ee*ee)+lam*np.trace(ee)**2/2)
    pc=np.atleast_1d(np.genfromtxt(p/'observables.csv',delimiter=',',names=True))
    rc=np.atleast_1d(np.genfromtxt(r/'observables.csv',delimiter=',',names=True))
    pl=pc[cycle*ps['steps_per_cycle']-1];rl=rc[cycle*rs['spc']-1]
    assert abs(pl['t']-rl['t'])<1e-10
    return dict(comparison_version=2,production=str(p.relative_to(ROOT)),reference=str(r.relative_to(ROOT)),cycle=cycle,
        production_dt=ps['dt'],reference_dt=8/rs['spc'],
        reference_lag_history=rs.get('lag_history',False),
        primary_L2_errors=dict(zip(('ux','uy','c','d'),np.sqrt(error).tolist())),
        primary_relative_L2_errors=dict(zip(('ux','uy','c','d'),np.sqrt(error/norm).tolist())),
        observable_relative_errors={k:float((pl[k]-rl[k])/rl[k]) for k in ('d_max','D_int','C_int','E_el')},
        host_Q1_energy=float(pl['E_el']),independently_recovered_enhanced_energy=float(enhanced_energy),
        reference_sampled_on_production_nodes_Q1_energy=reference_Q1_energy,
        reference_energy=float(rl['E_el']))


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--production');ap.add_argument('--reference');ap.add_argument('--cycle',type=int,default=5)
    ap.add_argument('--out',default='verification/field_comparison.json')
    args=ap.parse_args()
    result=dict(interpolation=check())
    if args.production or args.reference:
        assert args.production and args.reference
        result['comparison']=compare((ROOT/args.production).resolve(),(ROOT/args.reference).resolve(),args.cycle)
    (ROOT/args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
