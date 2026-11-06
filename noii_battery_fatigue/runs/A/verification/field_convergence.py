"""Compare actual Q1 fields at common physical points, including warped quads."""
import json
import numpy as np
from scipy.spatial import cKDTree
from common import ROOT,meshes,write_json
from manufactured import geometry


def interpolate(nodes,elems,uv,points,with_gradient=False):
    X=nodes[elems]
    tree=cKDTree(X.mean(axis=1))
    nearest=tree.query(points,k=min(16,len(elems)))[1]
    if nearest.ndim==1:
        nearest=nearest[:,None]
    s=np.array([-1.,1.,1.,-1.]); t=np.array([-1.,-1.,1.,1.])
    found=np.zeros(len(points),bool)
    result=np.zeros((len(points),uv.shape[1]))
    gradient=np.zeros((len(points),uv.shape[1],2))
    for col in range(nearest.shape[1]):
        todo=np.flatnonzero(~found)
        if not len(todo):
            break
        ids=nearest[todo,col]; local=X[ids]
        parent=np.zeros((len(todo),2))
        good=np.ones(len(todo),bool)
        for _ in range(12):
            xi,eta=parent.T
            N=.25*(1+s[None,:]*xi[:,None])*(1+t[None,:]*eta[:,None])
            dn=np.stack([.25*s[None,:]*(1+t[None,:]*eta[:,None]),
                         .25*t[None,:]*(1+s[None,:]*xi[:,None])],axis=-1)
            error=points[todo]-np.einsum('ma,mai->mi',N,local)
            J=np.einsum('mai,maj->mij',local,dn)
            det=J[:,0,0]*J[:,1,1]-J[:,0,1]*J[:,1,0]
            good &= np.abs(det)>1e-16
            det=np.where(np.abs(det)>1e-16,det,1.)
            update=np.column_stack([(J[:,1,1]*error[:,0]-J[:,0,1]*error[:,1])/det,
                                    (-J[:,1,0]*error[:,0]+J[:,0,0]*error[:,1])/det])
            parent=np.clip(parent+update,-5.,5.)
        xi,eta=parent.T
        N=.25*(1+s[None,:]*xi[:,None])*(1+t[None,:]*eta[:,None])
        error=points[todo]-np.einsum('ma,mai->mi',N,local)
        inside=good&(np.max(np.abs(parent),axis=1)<1+1e-10)&(np.linalg.norm(error,axis=1)<1e-11)
        result[todo[inside]]=np.einsum('ma,maj->mj',N[inside],uv[elems[ids[inside]]])
        if with_gradient:
            dn=np.stack([.25*s[None,:]*(1+t[None,:]*eta[:,None]),
                         .25*t[None,:]*(1+s[None,:]*xi[:,None])],axis=-1)
            J=np.einsum('mai,maj->mij',local[inside],dn[inside])
            B=np.einsum('mai,mij->maj',dn[inside],np.linalg.inv(J))
            gradient[todo[inside]]=np.einsum('mak,maj->mjk',B,uv[elems[ids[inside]]])
        found[todo[inside]]=True
    if not found.all():
        raise ValueError(f'Could not map {np.count_nonzero(~found)} quadrature points')
    return (result,gradient) if with_gradient else result


def samples(path):
    a=np.load(path); uv=a['U'].reshape(len(a['nodes']),-1)
    pts=[]; weights=[]; val=[]
    for N,B,w,points in geometry(a['nodes'],a['elems'],3):
        pts.append(points); weights.append(w)
        val.append(np.einsum('a,maj->mj',N,uv[a['elems']]))
    return a,uv,np.concatenate(pts),np.concatenate(weights),np.concatenate(val)


def norms(values,weights):
    squared=np.einsum('m,mj->j',weights,values**2)
    return {'u':float(np.sqrt(squared[:2].sum())),
            'c':float(np.sqrt(squared[2])),'d':float(np.sqrt(squared[3])),
            'q':float(np.sqrt(squared[4]))}


def main():
    # Check the inverse-map sampler against exactly reproduced affine fields.
    x,e=meshes.distorted(8,seed=4)
    uv=np.column_stack([2*x[:,0]+3*x[:,1]-.2,-x[:,0]+x[:,1]])
    points=np.random.default_rng(11).uniform(0,1,(1000,2))
    truth=np.column_stack([2*points[:,0]+3*points[:,1]-.2,-points[:,0]+points[:,1]])
    sampled,gradient=interpolate(x,e,uv,points,with_gradient=True)
    sampler_error=float(np.max(np.abs(sampled-truth)))
    gradient_error=float(np.max(np.abs(gradient-np.array([[2.,3.],[-1.,1.]])[None,:,:])))
    assert sampler_error<1e-12 and gradient_error<1e-12

    base=ROOT/'verification'/'runs'
    refpath=base/'selected_D0_s64_spc40_5cyc'/'fields_cycle_005.npz'
    ref,uv,pts,w,values=samples(refpath)
    spatial=[]
    for mesh in ['s','d']:
        for n in [8,16,32]:
            a=np.load(base/f'selected_D0_{mesh}{n}_spc40'/'fields_cycle_005.npz')
            field=a['U'].reshape(len(a['nodes']),6)
            val=interpolate(a['nodes'],a['elems'],field,pts)
            spatial.append({'mesh':'structured' if mesh=='s' else 'distorted','n':n,
                            'L2_difference_from_n64':norms(val-values,w)})
    temporal=[]
    for mesh in ['s','d']:
        for coarse,fine in [(40,80),(80,160)]:
            path=base/f'selected_D0_{mesh}16_spc{fine}'/'fields_cycle_020.npz'
            ref,uv,pts,w,val=samples(path)
            a=np.load(base/f'selected_D0_{mesh}16_spc{coarse}'/'fields_cycle_020.npz')
            other=interpolate(a['nodes'],a['elems'],a['U'].reshape(-1,6),pts)
            temporal.append({'mesh':'structured' if mesh=='s' else 'distorted',
                'coarse_spc':coarse,'fine_spc':fine,'L2_change':norms(other-val,w),
                'H_max_change':float(np.max(np.abs(ref['svars'][:,0::3]-a['svars'][:,0::3]))),
                'alpha_max_change':float(np.max(np.abs(ref['svars'][:,1::3]-a['svars'][:,1::3])))})
    result={'sampler_affine_max_error':sampler_error,'sampler_affine_gradient_error':gradient_error,'spatial_cycle_5':spatial,
            'temporal_cycle_20':temporal,
            'reference':'same selected UEL, structured n64 spc40, cycle5; not independent 2D reference'}
    write_json(ROOT/'verification'/'field_convergence.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
