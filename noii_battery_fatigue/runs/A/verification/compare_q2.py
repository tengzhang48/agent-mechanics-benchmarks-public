"""Compare the delivered Q1 method with the separate direct-Hessian Q2 family."""
import json
import numpy as np
from numpy.polynomial.legendre import leggauss
from common import ROOT,write_json
from field_convergence import interpolate
from q2_reference import shape


def q2_values(a,n,points):
    i=np.minimum((points[:,0]*n).astype(int),n-1)
    j=np.minimum((points[:,1]*n).astype(int),n-1)
    xi=2*(points[:,0]*n-i)-1;eta=2*(points[:,1]*n-j)-1
    lx=np.column_stack([xi*(xi-1)/2,1-xi**2,xi*(xi+1)/2])
    ly=np.column_stack([eta*(eta-1)/2,1-eta**2,eta*(eta+1)/2])
    N=(ly[:,:,None]*lx[:,None,:]).reshape(-1,9)
    uv=a['U'].reshape(-1,4)[a['elems'][j*n+i]]
    return np.einsum('ma,maj->mj',N,uv)


def norm(values,w):
    square=np.einsum('m,mj->j',w,values**2)
    return {'u':float(np.sqrt(square[:2].sum())),'c':float(np.sqrt(square[2])),
            'd':float(np.sqrt(square[3]))}


def main():
    out=ROOT/'verification';n=32
    ref=np.load(out/'q2_reference'/f'n{n}_spc40'/'cycle5.npz')
    uv=ref['U'].reshape(-1,4)[ref['elems']]
    gp,gw=leggauss(5);pts=[];weights=[];values=[];flux=[];qgrad=[]
    E=10.;nu=.28;mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu));pc=2*(mu+lam)/3
    beta=.035;c0=.35;kap=1e-5;D=.14;chi=.25
    for iy,eta in enumerate(gp):
        for ix,xi in enumerate(gp):
            N,B=shape(xi,eta,1/n)
            val=np.einsum('a,maj->mj',N,uv)
            grad=np.einsum('ak,maj->mjk',B,uv)
            lx=np.array([xi*(xi-1)/2,1-xi*xi,xi*(xi+1)/2]);ly=np.array([eta*(eta-1)/2,1-eta*eta,eta*(eta+1)/2])
            dx=np.array([xi-.5,-2*xi,xi+.5]);dy=np.array([eta-.5,-2*eta,eta+.5]);dd=np.array([1.,-2.,1.])
            hs=np.column_stack([np.outer(ly,dd).ravel(),np.outer(dy,dx).ravel(),np.outer(dd,lx).ravel()])*4*n*n
            hess=np.einsum('ak,maj->mjk',hs,uv[:,:,:2])
            pressure=pc*(grad[:,0,0]+grad[:,1,1]-2*beta*(val[:,2]-c0))
            pg=pc*(np.column_stack([hess[:,0,0]+hess[:,1,1],hess[:,0,1]+hess[:,1,2]])-2*beta*grad[:,2,:])
            g=(1-kap)*(1-val[:,3])**2+kap;gd=-2*(1-kap)*(1-val[:,3])
            qg=g[:,None]*pg+(pressure*gd)[:,None]*grad[:,3,:]
            J=(g*D)[:,None]*grad[:,2,:]-(chi*val[:,2])[:,None]*qg
            flux.append(J);qgrad.append(qg)
            pts.append(np.einsum('a,mai->mi',N,ref['nodes'][ref['elems']]))
            weights.append(np.full(len(uv),gw[ix]*gw[iy]/(4*n*n)))
            values.append(val)
    points=np.concatenate(pts);w=np.concatenate(weights);true=np.concatenate(values)
    flux=np.concatenate(flux);qgrad=np.concatenate(qgrad)
    q1=[];base=out/'runs'
    reference_history=np.load(out/'q2_reference'/f'n{n}_spc40'/'observables.npz')['history']
    for mesh,ns in [('s',[8,16,32,64]),('d',[8,16,32])]:
        for size in ns:
            suffix='_5cyc' if size==64 else ''
            path=base/f'selected_D0_{mesh}{size}_spc40{suffix}'
            a=np.load(path/'fields_cycle_005.npz')
            vals,grads=interpolate(a['nodes'],a['elems'],a['U'].reshape(-1,6),points,with_gradient=True)
            g=(1-kap)*(1-vals[:,3])**2+kap
            J=(D*g)[:,None]*grads[:,2,:]-(chi*vals[:,2])[:,None]*grads[:,4,:]
            hist=np.genfromtxt(path/'observables.csv',names=True,delimiter=',')[:200]
            obs=np.column_stack([hist[k] for k in ['d_max','D_int','C_int','E_el']])
            q1.append({'mesh':'structured' if mesh=='s' else 'distorted','n':size,
                'L2_field_difference':norm(vals[:,:4]-true,w),
                'flux_L2_difference':float(np.sqrt(np.sum(w[:,None]*(J-flux)**2))),
                'grad_q_L2_difference':float(np.sqrt(np.sum(w[:,None]*(grads[:,4,:]-qgrad)**2))),
                'max_abs_observable_difference':dict(zip(['d_max','D_int','C_int','E_el'],
                    np.max(np.abs(obs-reference_history[:,1:]),axis=0).tolist()))})
    selfcheck=[]
    for size in [4,8,16]:
        path=out/'q2_reference'/f'n{size}_spc40'
        a=np.load(path/'cycle5.npz');val=q2_values(a,size,points)
        hist=np.load(path/'observables.npz')['history']
        selfcheck.append({'n':size,'L2_field_difference':norm(val-true,w),
            'max_abs_observable_difference':dict(zip(['d_max','D_int','C_int','E_el'],
                np.max(np.abs(hist[:,1:]-reference_history[:,1:]),axis=0).tolist()))})
    result={'reference':'self-authored Q2 n32 direct Hessian, 5 cycles, spc40; same history time split',
            'reference_field_L2':norm(true,w),'Q1_comparison':q1,'Q2_self_convergence':selfcheck,
            'reference_flux_L2':float(np.sqrt(np.sum(w[:,None]*flux**2)))}
    write_json(out/'q2_comparison.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
