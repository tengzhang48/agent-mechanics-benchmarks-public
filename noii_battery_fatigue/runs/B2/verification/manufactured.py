"""Independent forced smooth solution: full mechanics/transport/Helmholtz PDE.

SymPy differentiates analytical primary fields to construct sources directly
from PUBLIC_TASK.md. The UEL sees no special code or fields. Verification-only
batch assembly subtracts these integrated loads; all solves use the host.
History H=0.7 is prescribed as an accepted prehistory and remains inactive.
"""
from pathlib import Path
import sys
import json
import argparse
import time
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.common import compiled, properties
from verification.reference import shape
import fehost
import meshes


def oracle(prm,dt,H0):
    x,y=sp.symbols("x y",real=True)
    u=sp.Matrix([.008*sp.sin(sp.pi*x)*sp.sin(sp.pi*y),
                 .006*sp.cos(sp.pi*x)*sp.sin(sp.pi*y)])
    c=prm["c0"]+.1*sp.cos(sp.pi*x)*sp.cos(sp.pi*y)
    d=.3+.1*sp.sin(sp.pi*x)*sp.sin(sp.pi*y)
    mu=prm["E"]/(2*(1+prm["nu"]))
    lam=prm["E"]*prm["nu"]/((1+prm["nu"])*(1-2*prm["nu"]))
    gu=u.jacobian([x,y]); ee=(gu+gu.T)/2-prm["beta"]*(c-prm["c0"])*sp.eye(2)
    g=(1-prm["kappa"])*(1-d)**2+prm["kappa"]
    stress=g*(2*mu*ee+lam*sp.trace(ee)*sp.eye(2))
    p=sp.trace(stress)/3
    grad=lambda f:sp.Matrix([sp.diff(f,x),sp.diff(f,y)])
    div=lambda a:sp.diff(a[0],x)+sp.diff(a[1],y)
    J=g*prm["D"]*grad(c)-prm["chi"]*c*grad(p)
    b=[-div(stress[i,:]) for i in range(2)]
    sc=(c-prm["c0"])/dt-div(J)
    sd=-prm["ell"]**2*div(grad(d))+d-(1-prm["kappa"])*(1-d)*H0
    fields=list(u)+[c,d,p,1-d]
    loads=b+[sc,-sd,sp.S(0),sp.S(0)]
    exact=sp.lambdify((x,y),fields,modules="numpy",cse=True)
    source=sp.lambdify((x,y),loads,modules="numpy",cse=True)
    gradients=sp.lambdify((x,y),[sp.diff(f,z) for f in fields for z in (x,y)],
                           modules="numpy",cse=True)
    def evaluate(fun,X):
        values=fun(X[:,0],X[:,1])
        return np.column_stack([np.broadcast_to(v,len(X)) for v in values])
    return tuple(lambda X,f=f:evaluate(f,X) for f in (exact,source,gradients))


def run(n,mesh,seed=1,contractfile=None,outdir=None):
    start=time.monotonic()
    contract=json.loads(Path(contractfile).read_text()) if contractfile else fehost.load_contract(str(ROOT))
    prm=properties(contract);prm.update(G_c=10.,alpha_T=1.)
    contract={**contract,"props":[prm[k] for k in contract["prop_names"]]}
    m=compiled(contract)
    nodes,elems=meshes.structured(n) if mesh=="structured" else meshes.distorted(n,seed)
    dt=.2;H0=.7
    exact,source,grad=oracle(prm,dt,H0)
    el=fehost.ContractUEL(m,contract,len(elems),dt)
    el.svars[:,0::3]=H0
    loads=np.zeros((len(elems),4,6))
    gp,gw=np.polynomial.legendre.leggauss(4)
    X=nodes[elems]
    for i,xi in enumerate(gp):
        for j,eta in enumerate(gp):
            N,dn=shape(xi,eta)
            J=np.einsum("mai,aj->mij",X,dn)
            weights=np.linalg.det(J)*gw[i]*gw[j]
            src=source(np.einsum("mai,a->mi",X,N))
            loads+=N[None,:,None]*src[:,None,:]*weights[:,None,None]
    loads=loads.reshape(len(elems),24)
    def batch(coords,U,DU):
        R,K=el.element_rk_batch(coords,U,DU)
        return R-loads,K
    bound=np.where((np.isclose(nodes,0)|np.isclose(nodes,1)).any(axis=1))[0]
    target=exact(nodes)
    bc={6*k+j:target[k,j] for k in bound for j in (0,1,2,3)}
    U0=np.zeros(len(nodes)*6);U0[2::6]=prm["c0"];U0[5::6]=1.
    U,report=fehost.solve(nodes,elems,el,lambda frac:bc,1,dt,U0=U0,batch_fn=batch,
        residual_scale=fehost.residual_weights(el,contract))
    values=U.reshape(-1,6)[elems]
    errors=np.zeros(6); gp_error=0.;flux_error=0.;norm_flux=0.
    gp,gw=np.polynomial.legendre.leggauss(5)
    for i,xi in enumerate(gp):
        for j,eta in enumerate(gp):
            N,dn=shape(xi,eta)
            J=np.einsum("mai,aj->mij",X,dn)
            B=np.einsum("aj,mjk->mak",dn,np.linalg.inv(J))
            weights=np.linalg.det(J)*gw[i]*gw[j]
            pts=np.einsum("mai,a->mi",X,N)
            ref=exact(pts)
            vf=np.einsum("maj,a->mj",values,N)
            gradients=np.einsum("maj,mak->mjk",values,B)
            rg=grad(pts).reshape(-1,6,2)
            errors+=np.sum((vf-ref)**2*weights[:,None],axis=0)
            gp_error+=np.sum((gradients[:,4]-rg[:,4])**2*weights[:,None])
            gf=(1-prm["kappa"])*vf[:,5]**2+prm["kappa"]
            gr=(1-prm["kappa"])*(1-ref[:,3])**2+prm["kappa"]
            jf=gf[:,None]*prm["D"]*gradients[:,2]-prm["chi"]*vf[:,2,None]*gradients[:,4]
            jr=gr[:,None]*prm["D"]*rg[:,2]-prm["chi"]*ref[:,2,None]*rg[:,4]
            flux_error+=np.sum((jf-jr)**2*weights[:,None])
            norm_flux+=np.sum(jr**2*weights[:,None])
    assert np.allclose(el.svars[:,0::3],H0,rtol=0,atol=0),"Manufactured H branch activated"
    dplus=U[3::6]+U[5::6]
    result=dict(n=n,mesh=mesh,seed=seed,wall_s=time.monotonic()-start,
        L2=dict(zip(contract["dof_names"],np.sqrt(errors).tolist())),
        grad_p_L2=float(np.sqrt(gp_error)),flux_L2=float(np.sqrt(flux_error)),
        flux_relative_L2=float(np.sqrt(flux_error/norm_flux)),
        complement_max_error=float(np.max(abs(dplus-1))),
        H_exactly_frozen=True,solver=report.summary())
    out=Path(outdir) if outdir else ROOT/"verification/manufactured"
    out.mkdir(parents=True,exist_ok=True)
    (out/f"{mesh}_{n}_s{seed}.json").write_text(json.dumps(result,indent=2)+"\n")
    np.savez_compressed(out/f"{mesh}_{n}_s{seed}.npz",nodes=nodes,elems=elems,U=U,
                        exact_nodal=target,dof_names=np.array(contract["dof_names"]))
    print(json.dumps(result),flush=True)
    return result


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--n",type=int,nargs="+",default=[8,16,32])
    ap.add_argument("--mesh",nargs="+",default=["structured","distorted"])
    ap.add_argument("--contract")
    ap.add_argument("--out")
    args=ap.parse_args()
    for mesh in args.mesh:
        for n in args.n:
            run(n,mesh,contractfile=args.contract,outdir=args.out)
