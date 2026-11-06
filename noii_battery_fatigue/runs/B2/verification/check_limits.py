"""Quantitative global oracles: Fourier diffusion and exhausted-fatigue affine patch."""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from verification.common import compiled, properties
from verification.reference import shape, constitutive
import fehost
import meshes


def mass_integral(nodes,elems,c):
    result=0.
    q,w=np.polynomial.legendre.leggauss(3)
    X=nodes[elems]
    for i,xi in enumerate(q):
        for j,eta in enumerate(q):
            N,dn=shape(xi,eta)
            det=np.linalg.det(np.einsum("mai,aj->mij",X,dn))
            result+=np.sum((c[elems]@N)*det)*w[i]*w[j]
    return float(result)


def diffusion(n,dt,mesh="structured"):
    c=fehost.load_contract(str(ROOT));prm=properties(c);prm["beta"]=0.
    c={**c,"props":[prm[p] for p in c["prop_names"]]}
    module=compiled(c)
    nodes,elems=meshes.structured(n) if mesh=="structured" else meshes.distorted(n,2)
    el=fehost.ContractUEL(module,c,len(elems),dt)
    steps=round(.5/dt);end=steps*dt
    mode=np.cos(np.pi*nodes[:,0])*np.cos(np.pi*nodes[:,1])
    U0=np.zeros(len(nodes)*6);U0[2::6]=prm["c0"]+.05*mode;U0[5::6]=1.
    n0=meshes.node_at(nodes,0,0);n1=meshes.node_at(nodes,1,0)
    bc={6*n0:0.,6*n0+1:0.,6*n1+1:0.}
    m0=mass_integral(nodes,elems,U0[2::6])
    changes=[]
    def output(step,frac,U,metrics=None):
        changes.append(mass_integral(nodes,elems,U[2::6])-m0)
    U,report=fehost.solve(nodes,elems,el,lambda frac:bc,steps,dt,U0=U0,
        residual_scale=fehost.residual_weights(el,c),step_callback=output)
    kh=6*(1-np.cos(np.pi/n))/((1/n)**2*(2+np.cos(np.pi/n)))
    afe=.05*(1+dt*prm["D"]*2*kh)**(-steps)
    abe=.05*(1+dt*prm["D"]*2*np.pi**2)**(-steps)
    acontinuous=.05*np.exp(-prm["D"]*2*np.pi**2*end)
    aobserved=float(mode@(U[2::6]-prm["c0"])/(mode@mode))
    error=np.max(abs(U[2::6]-(prm["c0"]+afe*mode)))
    if mesh=="structured":
        assert error<2e-11,error
    assert max(abs(np.array(changes)))<2e-12
    assert np.max(abs(U[3::6]))<1e-12
    assert np.max(abs(el.svars))<1e-12
    # Integrated spatial error against the exact backward-Euler continuum mode.
    q,w=np.polynomial.legendre.leggauss(4);l2=0.
    for i,xi in enumerate(q):
        for j,eta in enumerate(q):
            N,dn=shape(xi,eta);X=nodes[elems]
            det=np.linalg.det(np.einsum("mai,aj->mij",X,dn))
            pts=np.einsum("mai,a->mi",X,N)
            exact=prm["c0"]+abe*np.cos(np.pi*pts[:,0])*np.cos(np.pi*pts[:,1])
            err=U[2::6][elems]@N-exact
            l2+=np.sum(err**2*det)*w[i]*w[j]
    return dict(n=n,mesh=mesh,dt=dt,steps=steps,observed_amplitude=aobserved,
        FE_exact_amplitude=afe,BE_exact_amplitude=abe,continuous_amplitude=acontinuous,
        nodal_FE_error=float(error),L2_spatial_error=float(np.sqrt(l2)),
        max_mass_change=float(max(abs(np.array(changes)))),solver=report.summary())


def near_failure(n=16):
    c=fehost.load_contract(str(ROOT));prm=properties(c);module=compiled(c)
    nodes,elems=meshes.distorted(n,seed=5)
    el=fehost.ContractUEL(module,c,len(elems),1.)
    H0=1e12;el.svars[:,0::3]=H0;el.svars[:,1::3]=.02
    G=np.array([[.004,.001],[.002,.002]])
    target=nodes@G.T
    outer=np.where((np.isclose(nodes,0)|np.isclose(nodes,1)).any(axis=1))[0]
    bc={6*k+j:target[k,j] for k in outer for j in (0,1)}
    bc.update({6*k+2:prm["c0"] for k in range(len(nodes))})
    s0=1/(1+(1-prm["kappa"])*H0)
    U0=np.zeros((len(nodes),6));U0[:,2]=prm["c0"];U0[:,3]=1-s0;U0[:,5]=s0
    U,report=fehost.solve(nodes,elems,el,lambda frac:bc,2,.5,U0=U0.ravel(),
        residual_scale=fehost.residual_weights(el,c))
    # At fatigue floor, degradation differs from kappa by O(1/H^2).
    stress,pressure,state,g=constitutive(G,prm["c0"],1.,np.array([H0,.02,0.]),prm)
    H=state[0];sexact=1/(1+(1-prm["kappa"])*H)
    values=U.reshape(-1,6)
    ue=float(np.max(abs(values[:,:2]-target)))
    se=float(np.max(abs(values[:,5]-sexact))/sexact)
    he=float(np.max(abs(el.svars[:,0::3]-H))/H)
    pe=float(np.max(abs(values[:,4]-pressure)))
    assert ue<2e-10,(ue,se,he,pe)
    assert se<1e-7,(ue,se,he,pe)
    assert he<1e-7,(ue,se,he,pe)
    assert pe<1e-11,(ue,se,he,pe)
    return dict(n=n,mesh="distorted",seed=5,H_exact=float(H),s_exact=float(sexact),
        displacement_max_error=ue,s_max_relative_error=se,H_max_relative_error=he,
        pressure_max_error=pe,d_min=float(values[:,3].min()),
        d_max=float(values[:,3].max()),solver=report.summary())


def original_cancellation_control():
    # Retained rejected formulation, freshly called through the same host.
    c=json.loads((ROOT/"verification/exploration/initial5/uel_contract.json").read_text())
    c={**c,"source":["verification/exploration/initial5/generated_kernel.for",
                    "verification/exploration/initial5/uel.f90"]}
    m=compiled(c);el=fehost.ContractUEL(m,c,1,.1)
    X=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
    q=np.zeros((4,5));q[:,0]=.01*X[:,1];q[:,2]=.35
    Uold=q.copy().ravel();el.svars[:]=np.tile([0.,.02,0.],4)
    trace=[]
    for iteration in range(8):
        R,K=el.element_rk(0,X,q.ravel(),q.ravel()-Uold)
        trace.append(float(np.linalg.norm(R[3::5])))
        q[:,3]+=np.linalg.solve(K[3::5,3::5],-R[3::5])
    assert min(trace[1:])>1e-3
    return dict(residual_norms=trace,status="known-broken cancellation detected")


def main():
    result={"status":"running"}
    result["diffusion_spatial"]=[diffusion(n,.05,mesh) for mesh in ("structured","distorted") for n in (8,16,32)]
    result["diffusion_time"]=[diffusion(16,dt) for dt in (.1,.05,.025,.0125)]
    result["fully_damaged_global_patch"]=[near_failure(n) for n in (8,16,32)]
    result["old_damage_cancellation_control"]=original_cancellation_control()
    result["status"]="passed"
    (ROOT/"verification/limit_checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
