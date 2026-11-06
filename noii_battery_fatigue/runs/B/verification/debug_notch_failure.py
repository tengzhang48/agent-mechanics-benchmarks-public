"""Read-only production diagnosis of a captured stalled notch increment.

Experimental sources are copied/compiled under workspace tmp only. The main
element, generator, and host are never changed by this script.
"""
import importlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigs, spsolve

from common import ROOT, compiled, fehost
from check_element import quadrature


def main():
    capture = ROOT/"verification/custom_runs/notch_distorted_16_debug/failed_trial.npz"
    z = np.load(capture)
    nodes, elems = z["nodes"], z["elems"]
    ue, due, old = z["Ue"], z["DUe"], z["svars"]
    m, contract = compiled()
    parameters = json.loads((ROOT/"verification/problems/notch_sides.json").read_text())["parameters"]
    contract = {**contract, "props": [parameters.get(k,v) for k,v in zip(contract["prop_names"],contract["props"])]}
    ndof, dpn = len(nodes)*6, 6
    gm = (elems[:,:,None]*dpn+np.arange(dpn)).reshape(len(elems),24)
    x,y = nodes.T
    constrained = np.zeros(ndof,bool)
    constrained[6*np.where(np.isclose(x,0)|np.isclose(x,1))[0]+2]=True
    constrained[6*np.where(np.isclose(y,.5)&(x<=.125+1e-12))[0]+3]=True
    for xx in (0.,1.):
        k=np.where(np.isclose(x,xx)&np.isclose(y,0))[0][0]
        constrained[[6*k,6*k+1]]=True
    free=~constrained
    weights=np.tile([100.,100.,1.,1.,1000.,1000.],len(nodes))

    def assemble(rr,kk):
        residual=np.zeros(ndof)
        np.add.at(residual,gm.ravel(),rr.ravel())
        row=np.broadcast_to(gm[:,:,None],kk.shape).ravel()
        col=np.broadcast_to(gm[:,None,:],kk.shape).ravel()
        tangent=coo_matrix((kk.ravel(),(row,col)),shape=(ndof,ndof)).tocsr()
        return residual,tangent

    def element(module):
        el=fehost.ContractUEL(module,contract,len(elems),float(z["dt"]))
        el._time=z["time"]
        el.svars=old.copy()
        return el

    base_el=element(m)
    rr,kk=base_el.element_rk_batch(nodes[elems],ue,due)
    assert np.array_equal(rr,z["Re"]) and np.array_equal(kk,z["Ke"])
    R,K=assemble(rr,kk)
    du=np.zeros(ndof)
    du[free]=spsolve(K[free][:,free],-R[free])
    Rnorm=float(np.linalg.norm((weights*R)[free]))

    def norm_at(el,step):
        re,ke=el.element_rk_batch(nodes[elems],ue+step[gm],due+step[gm])
        residual,_=assemble(re,ke)
        return float(np.linalg.norm((weights*residual)[free])),residual

    norm_full,R_full=norm_at(base_el,du)
    physical_direction_error=(R_full-R-K@du)
    Kw=diags(weights[free])@K[free][:,free]
    eigenvalues=eigs(Kw,k=6,sigma=0,return_eigenvectors=False)
    report=dict(time=float(z["time"][1]),dt=float(z["dt"]),
                residual=Rnorm,
                component_residual_norm=[float(np.linalg.norm(R[j::6][free[j::6]])) for j in range(6)],
                max_newton_step_by_field=[float(np.max(abs(du[j::6]))) for j in range(6)],
                linear_relative_residual=float(np.linalg.norm((K@du+R)[free])/np.linalg.norm(R[free])),
                smallest_weighted_eigenvalues=[[float(v.real),float(v.imag)] for v in eigenvalues],
                full_step_actual_residual=norm_full,
                full_step_linearization_error=float(np.linalg.norm(physical_direction_error[free])),
                baseline_line_search={str(a):norm_at(base_el,a*du)[0] for a in (1.,.5,.1,.01,.001,1e-5,1e-7)})

    # An independent history-direction calculation at every Gauss point.
    p=parameters; mu=p["E"]/(2*(1+p["nu"])); lam=p["E"]*p["nu"]/((1+p["nu"])*(1-2*p["nu"]))
    bulk=lam+2*mu/3; a=1-p["kappa"]
    predicted=np.zeros(ndof); point_info=[]
    for ei,cell in enumerate(elems):
        values=ue[ei].reshape(4,6); direction=du[gm[ei]].reshape(4,6)
        for gp,(N,B,w) in enumerate(quadrature(nodes[cell])):
            strain=(values[:,:2].T@B+B.T@values[:,:2])/2
            delta_strain=(direction[:,:2].T@B+B.T@direction[:,:2])/2
            c,d,q,v=N@values[:,2:]; dc,dd,dq,dv=N@direction[:,2:]
            ee=strain-p["beta"]*(c-p["c0"])*np.eye(2)
            de=delta_strain-p["beta"]*dc*np.eye(2)
            tr=np.trace(ee); dt=np.trace(de)
            psi=bulk/2*max(tr,0)**2+mu*((ee[0,0]-ee[1,1])**2/2+2*ee[0,1]**2)
            delta_psi=bulk*max(tr,0)*dt+mu*((ee[0,0]-ee[1,1])*(de[0,0]-de[1,1])+4*ee[0,1]*de[0,1])
            h0,alpha0,w0=old[ei].reshape(4,3)[gp]
            degradation=p["kappa"]+a*v*v; W=degradation*psi
            alpha=alpha0+max(0.,W-w0)
            f=1. if alpha<=p["alpha_T"] else max(1e-15,max(0.,1-p["k_fg"]*math.log10(alpha/p["alpha_T"]))**2)
            if f!=1.:
                raise AssertionError("this diagnostic expects the captured undiminished-fatigue case")
            candidate=2*p["ell"]*psi/(p["G_c"]*f)
            delta_candidate=2*p["ell"]*delta_psi/p["G_c"]
            gap=candidate-h0
            tolerance=64*np.finfo(float).eps*max(1.,abs(candidate),abs(h0))
            near=abs(gap)<=tolerance
            if near and gap<=0 and delta_candidate>0:
                contribution=-a*v*delta_candidate*N*w
                np.add.at(predicted,cell*6+3,contribution)
            if near:
                point_info.append(dict(element=ei,gp=gp,H_old=float(h0),candidate=float(candidate),
                                       gap=float(gap),delta_candidate=float(delta_candidate),
                                       intact=float(v),alpha=float(alpha)))
    report["near_history_ties"]=point_info
    report["tie_direction_prediction_error"]=float(np.linalg.norm((predicted-physical_direction_error)[free]))
    report["predicted_tie_direction_error_norm"]=float(np.linalg.norm(predicted[free]))
    def delta_candidate_at_problem(step):
        ei,gp=127,1
        N,B,_=list(quadrature(nodes[elems[ei]]))[gp]
        values=ue[ei].reshape(4,6); direction=step[gm[ei]].reshape(4,6)
        strain=(values[:,:2].T@B+B.T@values[:,:2])/2
        delta_strain=(direction[:,:2].T@B+B.T@direction[:,:2])/2
        ee=strain-p["beta"]*(N@values[:,2]-p["c0"])*np.eye(2)
        de=delta_strain-p["beta"]*(N@direction[:,2])*np.eye(2)
        delta_psi=bulk*max(np.trace(ee),0)*np.trace(de)+mu*((ee[0,0]-ee[1,1])*(de[0,0]-de[1,1])+4*ee[0,1]*de[0,1])
        return float(2*p["ell"]*delta_psi/p["G_c"])
    report["frozen_branch_delta_candidate"]=delta_candidate_at_problem(du)

    # Production versus small-CS and alternate equality policies. No source is
    # changed outside tmp/history_tangent_<name>/.
    source=(ROOT/"element/generated_model.for").read_text()
    report["variants"]={}
    for name,policy in (("cs30","none"),("cs30_ge","geH"),
                        ("roundH_cs10","roundH"),("geW_cs10","geW"),
                        ("geBoth_cs10","geBoth")):
        folder=ROOT/"tmp"/f"history_tangent_{name}"
        folder.mkdir(parents=True,exist_ok=True)
        experimental=(source.replace("CS_H = 1.0d-10","CS_H = 1.0d-30")
                      if name.startswith("cs30") else source)
        if policy in ("geH","geBoth"):
            old_if="IF (DBLE(candidate) .GT. DBLE(hhist_old)) THEN"
            assert experimental.count(old_if)==1
            experimental=experimental.replace(old_if,"IF (DBLE(candidate) .GE. DBLE(hhist_old)) THEN")
        if policy in ("geW","geBoth"):
            old_if="IF (DBLE(w_new) .GT. DBLE(wprev_old)) THEN"
            assert experimental.count(old_if)==1
            experimental=experimental.replace(old_if,"IF (DBLE(w_new) .GE. DBLE(wprev_old)) THEN")
        if policy=="roundH":
            old_branch=("      IF (DBLE(candidate) .GT. DBLE(hhist_old)) THEN\n"
                        "      hhist_new = candidate\n      END IF")
            new_branch=("      IF (DBLE(candidate) .GT. DBLE(hhist_old)) THEN\n"
                        "      hhist_new = candidate\n"
                        "      ELSE IF (DBLE(hhist_old)-DBLE(candidate) .LE.\n"
                        "     & 512d0*EPSILON(1d0)*MAX(1d0,ABS(DBLE(hhist_old)),\n"
                        "     & ABS(DBLE(candidate)))) THEN\n"
                        "      hhist_new=hhist_old+DCMPLX(0d0,AIMAG(candidate))\n"
                        "      END IF")
            assert experimental.count(old_branch)==1
            experimental=experimental.replace(old_branch,new_branch)
        (folder/"generated_model.for").write_text(experimental)
        module_name=f"notch_history_{name}"
        work=folder/"build"
        if list(work.glob(module_name+"*.so")):
            sys.path.insert(0,str(work)); module=importlib.import_module(module_name)
        else:
            cc={**contract,"source":[str((folder/"generated_model.for").relative_to(ROOT)),"element/uel.f90"]}
            module=fehost.build_module(str(ROOT),cc,str(work),module_name)
        el=element(module)
        re,ke=el.element_rk_batch(nodes[elems],ue,due)
        rv,kv=assemble(re,ke)
        step=np.zeros(ndof);step[free]=spsolve(kv[free][:,free],-rv[free])
        norm_new,rnew=norm_at(el,step)
        report["variants"][name]=dict(
            residual_bitwise_equal=bool(np.array_equal(re,rr)),
            element_K_max_difference=float(np.max(abs(ke-kk))),
            directional_linearization_error=float(np.linalg.norm((rnew-rv-kv@step)[free])),
            full_step_residual=norm_new,
            delta_candidate_at_problem=delta_candidate_at_problem(step),
            step_max_by_field=[float(np.max(abs(step[j::6]))) for j in range(6)])
        if policy=="roundH":
            active_tangent=kv
    # The alternate H choice is a rank-one change. Its Sherman-Morrison
    # denominator locates a singular generalized Jacobian between branches.
    Kf=K[free][:,free]
    difference=(active_tangent[free][:,free]-Kf).toarray()
    rr=np.where(np.any(difference!=0,axis=1))[0]
    cc=np.where(np.any(difference!=0,axis=0))[0]
    left,singular,right=np.linalg.svd(difference[np.ix_(rr,cc)],full_matrices=False)
    u=np.zeros(free.sum());v=u.copy()
    u[rr]=left[:,0]*math.sqrt(singular[0])
    v[cc]=right[0]*math.sqrt(singular[0])
    coefficient=float(v@spsolve(Kf,u))
    report["rank_one_kink"]=dict(singular_values=singular.tolist(),
                                  sherman_morrison_coefficient=coefficient,
                                  singular_generalized_jacobian_fraction=-1/coefficient)
    (ROOT/"verification/notch_failure_review.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
