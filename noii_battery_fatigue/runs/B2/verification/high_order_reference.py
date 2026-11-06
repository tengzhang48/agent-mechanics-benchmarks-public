"""Independent Python Qp reference with analytical physical shape Hessians.

This is NOT the submitted UEL, and does not use abaqus_ufl. It integrates the
public laws directly with continuous Qp u,c,d, p>=2, and the broken *strong*
gradient of g*p_eff. That gradient is consistent under high-order refinement;
no projected pressure or intactness alias is used. The unchanged host provides
only generic assembly, Newton, line search and transactional stepping.
Scope: resolved smooth D0, away from the numerical fatigue-exhaustion limit.
"""
from pathlib import Path
import sys
import json
import argparse
import time
from dataclasses import asdict
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"host"))
import fehost
import meshes


def basis(order,r,s):
    pts=np.linspace(-1,1,order+1)
    pol=[]
    for k in range(order+1):
        f=np.polynomial.Polynomial.fromroots(np.delete(pts,k))
        pol.append(f/f(pts[k]))
    a=np.array([f(r) for f in pol]);b=np.array([f(s) for f in pol])
    da=np.array([f.deriv()(r) for f in pol]);db=np.array([f.deriv()(s) for f in pol])
    dda=np.array([f.deriv(2)(r) for f in pol]);ddb=np.array([f.deriv(2)(s) for f in pol])
    N=np.outer(b,a).ravel()
    B=np.stack([np.outer(b,da).ravel(),np.outer(db,a).ravel()],axis=-1)
    H=np.empty((len(N),2,2))
    H[:,0,0]=np.outer(b,dda).ravel();H[:,1,1]=np.outer(ddb,a).ravel()
    H[:,0,1]=H[:,1,0]=np.outer(db,da).ravel()
    return N,B,H


def mesh_qp(n,order,distorted=False,seed=1):
    nodes,_=meshes.structured(n*order)
    coarse,quads=meshes.distorted(n,seed) if distorted else meshes.structured(n)
    elems=[]; grid=n*order+1
    for j in range(n):
        for i in range(n):
            conn=np.array([(order*j+b)*grid+order*i+a for b in range(order+1) for a in range(order+1)])
            elems.append(conn)
            X=coarse[quads[j*n+i]]
            for b,s in enumerate(np.linspace(-1,1,order+1)):
                for a,r in enumerate(np.linspace(-1,1,order+1)):
                    N=np.array([(1-r)*(1-s),(1+r)*(1-s),(1+r)*(1+s),(1-r)*(1+s)])/4
                    nodes[conn[b*(order+1)+a]]=N@X
    return nodes,np.array(elems)


class Reference:
    dof_per_node=4

    def __init__(self,X,order,prm,dt,lag_history=False,qorder=None):
        self.prm=prm;self.dt=dt;self._time=np.array([dt,dt]);self._kstep=1;self._kinc=1
        self.lag_history=lag_history
        self.order=order;self.nn=(order+1)**2;self.ndof=4*self.nn
        qorder=qorder or order+1
        ne=len(X);ng=qorder**2
        self.svars=np.zeros((ne,ng*3))
        self.Akin=np.zeros((ne,ng,16,self.ndof))
        self.Tests=np.zeros((ne,ng,self.ndof,9))
        self.weights=np.zeros((ne,ng));self.N=np.zeros((ng,self.nn));self.points=np.zeros((ne,ng,2))
        q,w=np.polynomial.legendre.leggauss(qorder)
        k=0
        for j,s in enumerate(q):
            for i,r in enumerate(q):
                N,dN,hN=basis(order,r,s)
                J=np.einsum("mai,ar->mir",X,dN)
                inv=np.linalg.inv(J)
                B=np.einsum("ar,mri->mai",dN,inv)
                hX=np.einsum("mai,ars->mirs",X,hN)
                corrected=hN[None]-np.einsum("mai,mirs->mars",B,hX)
                H=np.einsum("mri,mars,msj->maij",inv,corrected,inv)
                self.N[k]=N;self.points[:,k]=np.einsum("a,mai->mi",N,X)
                self.weights[:,k]=np.linalg.det(J)*w[i]*w[j]
                for a in range(self.nn):
                    ux,uy,c,d=range(4*a,4*a+4)
                    self.Akin[:,k,0,ux]=B[:,a,0];self.Akin[:,k,1,ux]=B[:,a,1]
                    self.Akin[:,k,2,uy]=B[:,a,0];self.Akin[:,k,3,uy]=B[:,a,1]
                    for off,comp in ((ux,4),(uy,7)):
                        self.Akin[:,k,comp,off]=H[:,a,0,0]
                        self.Akin[:,k,comp+1,off]=H[:,a,0,1]
                        self.Akin[:,k,comp+2,off]=H[:,a,1,1]
                    self.Akin[:,k,10,c]=N[a];self.Akin[:,k,11,d]=N[a]
                    self.Akin[:,k,12:14,c]=B[:,a,:]
                    self.Akin[:,k,14:16,d]=B[:,a,:]
                    self.Tests[:,k,ux,0]=B[:,a,0];self.Tests[:,k,ux,2]=B[:,a,1]
                    self.Tests[:,k,uy,1]=B[:,a,1];self.Tests[:,k,uy,2]=B[:,a,0]
                    self.Tests[:,k,c,3]=N[a];self.Tests[:,k,c,4:6]=B[:,a,:]
                    self.Tests[:,k,d,6]=N[a];self.Tests[:,k,d,7:9]=B[:,a,:]
                k+=1
        self.geom_H=H

    def law(self,z,old,cold):
        p=self.prm
        mu=p["E"]/(2*(1+p["nu"]))
        lam=p["E"]*p["nu"]/((1+p["nu"])*(1-2*p["nu"]))
        c,d=z[...,10],z[...,11]
        x=z[...,0]-p["beta"]*(c-p["c0"])
        y=z[...,3]-p["beta"]*(c-p["c0"])
        xy=(z[...,1]+z[...,2])/2
        tr=x+y
        g=(1-p["kappa"])*(1-d)**2+p["kappa"]
        sigma=np.stack((g*(2*mu*x+lam*tr),g*(2*mu*y+lam*tr),g*2*mu*xy),axis=-1)
        pos=np.where(tr.real>0,tr,0.)
        psi=(lam+2*mu/3)/2*pos**2+mu*((x-y)**2/2+2*xy**2)
        W=g*psi
        alpha=old[...,1]+np.where(W.real>old[...,2],W-old[...,2],0.)
        asafe=np.where(alpha.real>p["alpha_T"],alpha,p["alpha_T"])
        inside=1-p["k_fg"]*np.log(asafe/p["alpha_T"])/np.log(10.)
        fatigue=np.where(inside.real>0,inside,0.)**2
        fatigue=np.where(fatigue.real>1e-15,fatigue,1e-15)
        candidate=2*p["ell"]*psi/(p["G_c"]*fatigue)
        H=np.where(candidate.real>old[...,0],candidate,old[...,0])
        pe=2*(lam+mu)/3*tr
        gradpe=2*(lam+mu)/3*(np.stack((z[...,4]+z[...,8],z[...,5]+z[...,9]),axis=-1)
                                 -2*p["beta"]*z[...,12:14])
        gradg=-2*(1-p["kappa"])*(1-d[...,None])*z[...,14:16]
        gradp=g[...,None]*gradpe+pe[...,None]*gradg
        J=g[...,None]*p["D"]*z[...,12:14]-p["chi"]*c[...,None]*gradp
        drive=old[...,0] if self.lag_history else H
        features=np.concatenate((sigma,((c-cold)/self.dt)[...,None],J,
            (d-(1-p["kappa"])*(1-d)*drive)[...,None],p["ell"]**2*z[...,14:16]),axis=-1)
        state=np.stack((H,alpha,W),axis=-1)
        return features,state

    def element_rk_batch(self,coords,U,DU):
        z=np.einsum("mgza,ma->mgz",self.Akin,U,optimize=True)
        cold=np.einsum("ga,ma->mg",self.N,(U-DU)[:,2::4])
        old=self.svars.reshape(len(U),-1,3)
        features,state=self.law(z,old,cold)
        R=np.einsum("mgaf,mgf,mg->ma",self.Tests,features,self.weights,optimize=True)
        deriv=np.empty(z.shape[:2]+(9,16))
        for j in range(16):
            zz=z.astype(complex);zz[...,j]+=1e-30j
            deriv[...,j]=self.law(zz,old,cold)[0].imag/1e-30
        dFdu=np.einsum("mgfz,mgzb->mgfb",deriv,self.Akin,optimize=True)
        K=np.einsum("mgaf,mgfb,mg->mab",self.Tests,dFdu,self.weights,optimize=True)
        return R,K

    def element_rk(self,*args):
        raise AssertionError("Reference uses batched assembly only")

    def commit_step(self,nodes,elems,dpn,U,Uold):
        gm=(elems[:,:,None]*4+np.arange(4)).reshape(len(elems),-1)
        z=np.einsum("mgza,ma->mgz",self.Akin,U[gm],optimize=True)
        cold=np.einsum("ga,ma->mg",self.N,Uold[gm][:,2::4])
        new=self.law(z,self.svars.reshape(len(elems),-1,3),cold)[1]
        assert np.all(new[...,0]>=self.svars[:,0::3])
        assert np.all(new[...,1]>=self.svars[:,1::3])
        self.svars=new.reshape(len(elems),-1)

    def observables(self,nodes,elems,U):
        gm=(elems[:,:,None]*4+np.arange(4)).reshape(len(elems),-1)
        z=np.einsum("mgza,ma->mgz",self.Akin,U[gm],optimize=True)
        F,_=self.law(z,self.svars.reshape(len(elems),-1,3),z[...,10])
        p=self.prm
        xx=z[...,0]-p["beta"]*(z[...,10]-p["c0"])
        yy=z[...,3]-p["beta"]*(z[...,10]-p["c0"])
        xy=(z[...,1]+z[...,2])/2
        energy=(F[...,0]*xx+F[...,1]*yy+2*F[...,2]*xy)/2
        return dict(d_max=float(U[3::4].max()),D_int=float(np.sum(z[...,11]*self.weights)),
            C_int=float(np.sum(z[...,10]*self.weights)),E_el=float(np.sum(energy*self.weights)))


def check(qorder=None):
    prm=json.loads((ROOT/"host/dev_problem/d0.json").read_text())["parameters"]
    X=np.array([[0.,0.],[1.2,.1],[1.05,1.],[0.,1.1]])
    nodes=[]
    for y in np.linspace(-1,1,3):
        for x in np.linspace(-1,1,3):
            N=np.array([(1-x)*(1-y),(1+x)*(1-y),(1+x)*(1+y),(1-x)*(1+y)])/4
            nodes.append(N@X)
    nodes=np.array(nodes);ref=Reference(nodes[None],2,prm,.2,qorder=qorder)
    x,y=nodes.T
    U=np.column_stack((.01*x*x+.003*x*y,-.004*y*y+.002*x*y,.35+.1*x-.07*y,.2+.05*x+.03*y)).reshape(1,-1)
    z=np.einsum("mgza,ma->mgz",ref.Akin,U)
    pts=ref.points;x,y=pts[...,0],pts[...,1]
    expected=np.stack((.02*x+.003*y,.003*x,.002*y,-.008*y+.002*x,
        np.full_like(x,.02),np.full_like(x,.003),np.zeros_like(x),np.zeros_like(x),np.full_like(x,.002),np.full_like(x,-.008),
        .35+.1*x-.07*y,.2+.05*x+.03*y,np.full_like(x,.1),np.full_like(x,-.07),np.full_like(x,.05),np.full_like(x,.03)),axis=-1)
    geomerr=float(np.max(abs(z-expected)));assert geomerr<2e-14,geomerr
    area=abs(np.dot(X[:,0],np.roll(X[:,1],-1))-np.dot(X[:,1],np.roll(X[:,0],-1)))/2
    area_error=float(abs(ref.weights.sum()-area));assert area_error<2e-13,area_error
    ng=ref.weights.shape[1]
    old=np.zeros((1,ng,3));cold=np.full((1,ng),.35)
    feat,state=ref.law(z,old,cold)
    mu=prm["E"]/(2*(1+prm["nu"]));lam=prm["E"]*prm["nu"]/((1+prm["nu"])*(1-2*prm["nu"]))
    cc=expected[...,10];dd=expected[...,11];a=1-prm["kappa"]
    g=a*(1-dd)**2+prm["kappa"]
    pe=2*(mu+lam)/3*(.022*x-.005*y-2*prm["beta"]*(cc-.35))
    gradpe=2*(mu+lam)/3*(np.array([.022,-.005])-2*prm["beta"]*np.array([.1,-.07]))
    gradg=-2*a*(1-dd[...,None])*np.array([.05,.03])
    J=g[...,None]*prm["D"]*np.array([.1,-.07])-prm["chi"]*cc[...,None]*(g[...,None]*gradpe+pe[...,None]*gradg)
    je=float(np.max(abs(J-feat[...,4:6])));assert je<2e-14,je
    prev=np.zeros_like(U);prev[:,2::4]=.35
    R,K=ref.element_rk_batch(nodes[None],U,U-prev)
    eps=1e-8;fd=np.zeros_like(K)
    for j in range(U.shape[1]):
        step=np.zeros_like(U);step[:,j]=eps
        fd[:,:,j]=(ref.element_rk_batch(nodes[None],U+step,U+step-prev)[0]-
                    ref.element_rk_batch(nodes[None],U-step,U-step-prev)[0])/(2*eps)
    ke=float(np.max(abs(K-fd))/np.max(abs(K)));assert ke<2e-6,ke
    result=dict(status="passed",quadrature_order=qorder or 3,geometry_area_error=area_error,quadratic_physical_kinematics_max_error=geomerr,
        analytical_pressure_gradient_flux_error=je,assembled_tangent_relative_error=ke)
    filename=f"reference_checks_q{qorder}.json" if qorder else "reference_checks.json"
    (ROOT/"verification"/filename).write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result),flush=True)


def run(n=8,order=2,cycles=2,spc=40,distorted=False,problem=None,lag_history=False,qorder=None):
    start=time.monotonic();prob=json.loads(Path(problem or ROOT/"host/dev_problem/d0.json").read_text())
    assert prob["c_edges"]==["top"] and prob.get("notch") is None,"Reference scope is top charging without a notch"
    prm=prob["parameters"];dt=prob["loading"]["T"]/spc
    nodes,elems=mesh_qp(n,order,distorted)
    ref=Reference(nodes[elems],order,prm,dt,lag_history,qorder)
    top=np.where(np.isclose(nodes[:,1],1))[0]
    n0=meshes.node_at(nodes,0,0);n1=meshes.node_at(nodes,1,0)
    def bc(frac):
        cb=prm["c0"]+prob["loading"]["A"]*(1-np.cos(2*np.pi*frac*cycles))/2
        vals={4*k+2:cb for k in top};vals.update({4*n0:0.,4*n0+1:0.,4*n1+1:0.})
        return vals
    U0=np.zeros(len(nodes)*4);U0[2::4]=prm["c0"]
    name=f"{prob['id']}_Q{order}_{'distorted' if distorted else 'structured'}_{n}_spc{spc}_cy{cycles}"
    if lag_history:name+="_acceptedH"
    if qorder:name+=f"_q{qorder}"
    out=ROOT/"verification/reference_runs"/name;out.mkdir(parents=True,exist_ok=True)
    rows=[]
    def cb(step,frac,U,metrics=None):
        rows.append(dict(step=step,t=frac*cycles*prob["loading"]["T"],**ref.observables(nodes,elems,U)))
        if step%spc==0:
            np.savez_compressed(out/f"fields_cycle_{step//spc:03d}.npz",nodes=nodes,elems=elems,
                U=U,svars=ref.svars,dof_names=np.array(["ux","uy","c","d"]))
            print(name,"cycle",step//spc,rows[-1],flush=True)
    status="ok";report=None
    try:
        U,report=fehost.solve(nodes,elems,ref,bc,cycles*spc,dt,U0=U0,
            residual_scale={0:100.,1:100.,2:1.,3:1.},step_callback=cb)
    except Exception as exc:
        status=f"failed: {type(exc).__name__}: {exc}"
    keys=["step","t","d_max","D_int","C_int","E_el"]
    with (out/"observables.csv").open("w") as f:
        f.write(",".join(keys)+"\n")
        for row in rows:
            f.write(",".join(f"{row[k]:.16g}" for k in keys)+"\n")
    result=dict(status=status,n=n,order=order,cycles=cycles,spc=spc,distorted=distorted,
                wall_s=time.monotonic()-start,accepted_output_steps=len(rows),parameters=prm,
                lag_history=lag_history,quadrature_order=qorder or order+1)
    if report:
        result.update(report.summary())
        (out/"solver_report.json").write_text(json.dumps(asdict(report),indent=1)+"\n")
    (out/"run.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result),flush=True)


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--check",action="store_true")
    ap.add_argument("--n",type=int,default=8)
    ap.add_argument("--order",type=int,default=2)
    ap.add_argument("--cycles",type=int,default=2)
    ap.add_argument("--spc",type=int,default=40)
    ap.add_argument("--distorted",action="store_true")
    ap.add_argument("--problem")
    ap.add_argument("--lag-history",action="store_true")
    ap.add_argument("--qorder",type=int)
    args=vars(ap.parse_args())
    if args.pop("check"):
        check(qorder=args["qorder"])
    else:
        run(**args)
