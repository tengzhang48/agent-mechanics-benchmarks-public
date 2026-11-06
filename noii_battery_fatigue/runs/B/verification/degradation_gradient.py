"""Independent projection/flux refinement for the p_eff*grad(g) contribution.

Prescribe affine isotropic strain and nonuniform linear damage. Then pressure
is constant but the degraded pressure gradient is nonzero. Only q and vfrac
are free. This is an operator test with prescribed fields, not a free PDE run.
"""
import json

import numpy as np

from common import ROOT,compiled,fehost,source_identity
import meshes


def main():
    module,contract=compiled()
    prm=dict(zip(contract["prop_names"],contract["props"]))
    mu=prm["E"]/(2*(1+prm["nu"]))
    lam=prm["E"]*prm["nu"]/((1+prm["nu"])*(1-2*prm["nu"]))
    p=4*(mu+lam)*.02/3
    rows=[]
    for mesh in ("structured","distorted"):
        for n in (8,16,32):
            nodes,elems=meshes.structured(n) if mesh=="structured" else meshes.distorted(n,seed=4)
            el=fehost.ContractUEL(module,contract,len(elems),.1)
            dpn=el.dof_per_node
            initial=np.zeros((len(nodes),dpn));initial[:,el.comp("c")]=prm["c0"]
            bc={}
            for i,(x,y) in enumerate(nodes):
                for k,val in [(el.comp("ux"),.02*x),(el.comp("uy"),.02*y),
                              (el.comp("c"),prm["c0"]),(el.comp("d"),.2+.4*x)]:
                    bc[dpn*i+k]=val
            U,report=fehost.solve(nodes,elems,el,lambda _:bc,1,.1,U0=initial.ravel(),
                                   residual_scale=fehost.residual_weights(el,contract),tol=1e-10)
            uu=U.reshape(len(nodes),dpn)
            totals=np.zeros(4)
            pts,wts=np.polynomial.legendre.leggauss(4)
            # Explicit independent shape/Jacobian calculations, not generator helpers.
            s=np.array([-1.,1.,1.,-1.]);t=np.array([-1.,-1.,1.,1.])
            X=nodes[elems];qnod=uu[elems,el.comp("q")]
            for i,xi in enumerate(pts):
                for j,eta in enumerate(pts):
                    N=(1+s*xi)*(1+t*eta)/4
                    DN=np.column_stack([s*(1+t*eta),t*(1+s*xi)])/4
                    J=np.einsum("maj,ai->mji",X,DN)
                    B=np.einsum("ai,mij->maj",DN,np.linalg.inv(J))
                    x=X[:,:,0]@N
                    v=.8-.4*x
                    q_exact=p*(prm["kappa"]+(1-prm["kappa"])*v*v)
                    grad_exact=np.column_stack([-.8*p*(1-prm["kappa"])*v,np.zeros(len(x))])
                    qh=qnod@N;gh=np.einsum("ma,maj->mj",qnod,B)
                    w=np.linalg.det(J)*wts[i]*wts[j]
                    totals[0]+=np.sum((qh-q_exact)**2*w)
                    totals[1]+=np.sum(np.sum((gh-grad_exact)**2,axis=1)*w)
                    totals[2]+=np.sum(np.sum(grad_exact**2,axis=1)*w)
                    totals[3]+=np.sum(np.sum(gh**2,axis=1)*w)
            errs=np.sqrt(totals)
            assert errs[1]<errs[2]/4 and errs[3]>errs[2]/2
            constraint=np.max(abs((uu[:,el.comp("d")]-1)+uu[:,el.comp("vfrac")]))
            assert constraint<1e-12 and report.cutbacks==0
            rows.append(dict(mesh=mesh,n=n,q_L2=float(errs[0]),grad_q_L2=float(errs[1]),
                             flux_L2=float(prm["chi"]*prm["c0"]*errs[1]),
                             missing_product_gradient_error=float(errs[2]),
                             constraint_max=float(constraint),solver=report.summary()))
    for mesh in ("structured","distorted"):
        rr=[r for r in rows if r["mesh"]==mesh]
        assert all(a["q_L2"]>b["q_L2"] and a["grad_q_L2"]>b["grad_q_L2"] for a,b in zip(rr,rr[1:]))
    data=dict(scope="constant p_eff, nonconstant degradation; prescribed u,c,d",pressure=p,
              source_sha256=source_identity(),results=rows)
    (ROOT/"verification/degradation_gradient.json").write_text(json.dumps(data,indent=2)+"\n")
    print(json.dumps(data,indent=2))


if __name__=="__main__":
    main()
