"""Audit host outputs and quantify self-convergence without an answer key.

This reports differences to the finest tested discretization, not errors to
an independent solution of the fully coupled PDE. Independent exact oracles
are provided separately. No plotting packages are required or installed.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from common import ROOT

QUANTITIES=("d_max","D_int","C_int","E_el","A_damage")


def write_csv(path, rows):
    if not rows:
        # An interim cohort may not yet have its reference. Do not leave a
        # table from a previous source/version looking like a current result.
        path.unlink(missing_ok=True)
        return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)


def point_shape(xi,eta):
    s=np.array([-1.,1.,1.,-1.]);t=np.array([-1.,-1.,1.,1.])
    return ((1+s*xi)*(1+t*eta)/4,
            np.column_stack((s*(1+t*eta),t*(1+s*xi)))/4)


def independent_observables(data,params):
    """Independent matrix-strain implementation of the documented bridge."""
    nodes,elems,U=data["nodes"],data["elems"],data["U"]
    names=data["dof_names"].tolist()
    uu=U.reshape(len(nodes),len(names))[elems]
    E,nu=params["E"],params["nu"]
    mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu))
    X=nodes[elems]
    quantities={k:0. for k in ("D_int","C_int","E_el","A_damage")}
    cidx,didx=names.index("c"),names.index("d")
    uxidx,uyidx=names.index("ux"),names.index("uy")
    for xi in [-1/np.sqrt(3),1/np.sqrt(3)]:
        for eta in [-1/np.sqrt(3),1/np.sqrt(3)]:
            N,DN=point_shape(xi,eta)
            J=np.einsum("maj,ai->mji",X,DN)
            weight=np.linalg.det(J)
            B=np.einsum("ai,mij->maj",DN,np.linalg.inv(J))
            c=uu[:,:,cidx]@N;d=uu[:,:,didx]@N
            un=uu[:,:,[uxidx,uyidx]]
            gu=np.einsum("mai,maj->mij",un,B)
            ee=(gu+gu.transpose(0,2,1))/2
            ee-=params["beta"]*(c-params["c0"])[:,None,None]*np.eye(2)
            psi=mu*np.sum(ee*ee,axis=(1,2))+lam*np.trace(ee,axis1=1,axis2=2)**2/2
            g=params["kappa"]+(1-params["kappa"])*(1-d)**2
            quantities["D_int"]+=float(d@weight)
            quantities["C_int"]+=float(c@weight)
            quantities["E_el"]+=float((g*psi)@weight)
    for xi in np.arange(4)/2-.75:
        for eta in np.arange(4)/2-.75:
            N,DN=point_shape(xi,eta)
            w=np.linalg.det(np.einsum("maj,ai->mji",X,DN))/4
            quantities["A_damage"]+=float((uu[:,:,didx]@N>=.5)@w)
    quantities["d_max"]=float(np.max(U[didx::len(names)]))
    return quantities


def compare_observables(data,reference):
    out={}
    for k in QUANTITIES:
        rr=np.interp(data["t"],reference["t"],reference[k])
        error=float(np.max(abs(data[k]-rr)))
        out[k+"_max_abs_diff"]=error
        reference_peak=float(np.max(abs(rr)))
        out[k+"_max_rel_diff"]=error/reference_peak if reference_peak>0. else None
    return out


def fine_q1(xy,reference):
    """Evaluate a structured finest-grid Q1 field at arbitrary physical points."""
    nf=round(np.sqrt(len(reference["nodes"])))-1
    vals=reference["U"].reshape(nf+1,nf+1,-1)
    p=np.clip(xy*nf,0.,float(nf))
    ij=np.minimum(np.floor(p).astype(int),nf-1)
    rx,ry=(p-ij).T;i,j=ij.T
    return ((1-rx)[:,None]*(1-ry)[:,None]*vals[j,i]
           +rx[:,None]*(1-ry)[:,None]*vals[j,i+1]
           +rx[:,None]*ry[:,None]*vals[j+1,i+1]
           +(1-rx)[:,None]*ry[:,None]*vals[j+1,i])


def compare_fields(data,reference):
    names=data["dof_names"].tolist()
    assert names==reference["dof_names"].tolist()
    X=data["nodes"][data["elems"]]
    uu=data["U"].reshape(len(data["nodes"]),len(names))[data["elems"]]
    err=np.zeros(len(names));den=np.zeros(len(names))
    pts,wts=np.polynomial.legendre.leggauss(3)
    for i,xi in enumerate(pts):
        for j,eta in enumerate(pts):
            N,DN=point_shape(xi,eta)
            xy=np.einsum("a,maj->mj",N,X)
            w=np.linalg.det(np.einsum("maj,ai->mji",X,DN))*wts[i]*wts[j]
            val=np.einsum("a,mak->mk",N,uu)
            ref=fine_q1(xy,reference)
            err+=np.sum((val-ref)**2*w[:,None],axis=0)
            den+=np.sum(ref**2*w[:,None],axis=0)
    return {k+"_L2_abs":float(np.sqrt(e)) for k,e in zip(names,err)} | {
        k+"_L2_rel":float(np.sqrt(e/max(d,1e-30))) for k,e,d in zip(names,err,den)}


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--require-complete",action="store_true")
    ap.add_argument("--version",default="selected")
    a=ap.parse_args()
    contract=json.loads((ROOT/'element/uel_contract.json').read_text())
    dpn=contract['dof_per_node'];nsv=contract['nsvars']
    hcols=np.asarray(contract['history_slots']['H'])-1
    acols=np.asarray(contract['history_slots']['alpha'])-1
    state_stride=nsv//len(hcols)
    records=[];table=[];audits=[];obs={};states={}
    prob=json.loads((ROOT/"host/dev_problem/d0.json").read_text())
    for path in sorted((ROOT/"verification/runs"/a.version).glob("*/run.json")):
        rec=json.loads(path.read_text());label=path.parent.name;records.append((label,rec))
        row=dict(label=label,mesh=rec["mesh"],n=rec["n"],spc=rec["steps_per_cycle"],
                 cycles=rec["cycles"],status=rec["status"],cutbacks=rec.get("cutbacks"),
                 min_dt=rec.get("min_dt"),wall_s=rec["wall_s"])
        for name,digest in rec.get("element_source_sha256",{}).items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,(label,name,"changed source")
        if rec["accepted_output_steps"]:
            oo=np.atleast_1d(np.genfromtxt(path.parent/"observables.csv",delimiter=",",names=True))
            assert np.all(np.diff(oo["t"])>0) and np.isfinite(oo.tolist()).all()
            obs[label]=oo
            row.update({k+"_final":float(oo[k][-1]) for k in QUANTITIES})
            row["d_max_peak"]=float(np.max(oo["d_max"]))
        files=sorted(path.parent.glob("fields_cycle_*.npz"))
        if rec["status"]=="ok":
            assert len(files)==rec["cycles"]
            assert rec["accepted_output_steps"]==rec["cycles"]*rec["steps_per_cycle"]
            assert all(v is not None and v>=0 for v in rec["history_min_change"].values())
        for file in files:
            data=dict(np.load(file));nn,ne=len(data["nodes"]),len(data["elems"])
            assert nn==(rec["n"]+1)**2 and ne==rec["n"]**2
            assert data["U"].size==dpn*nn and data["svars"].shape==(ne,nsv)
            assert data["elems"].shape==(ne,4)
            assert data["elems"].min()>=0 and data["elems"].max()<nn
            assert np.isfinite(data["U"]).all() and np.isfinite(data["svars"]).all()
            assert len(set(data["dof_names"].tolist()))==dpn
            cycle=int(file.stem.split("_")[-1])
            computed=independent_observables(data,prob["parameters"])
            stored=obs[label][cycle*rec["steps_per_cycle"]-1]
            error=max(abs(computed[k]-stored[k])/max(abs(computed[k]),1e-12) for k in QUANTITIES)
            assert error<2e-9,(label,cycle,error)
            uv=data["U"].reshape(nn,dpn)
            constraint=float(np.max(abs((uv[:,3]-1)+uv[:,5])))
            gdiff=float(np.max(abs(uv[:,5]**2-(1-uv[:,3])**2)))
            audits.append(dict(label=label,cycle=cycle,expected_nodes=nn,observed_nodes=nn,
                               expected_elements=ne,observed_elements=ne,missing=0,nonfinite=0,
                               bridge_max_relative_error=error,intact_constraint_max=constraint,
                               degradation_difference_without_kappa_factor=gdiff,
                               d_min=float(uv[:,3].min()),d_max=float(uv[:,3].max())))
            states[(label,cycle)]=data
            # H and alpha are irreversible, so their final-cycle maxima also
            # cover prior accepted increments at these same stored points.
            hmax=float(np.max(data['svars'][:,hcols]))
            amax=float(np.max(data['svars'][:,acols]))
            p=prob['parameters']
            fmin=1. if amax<=p['alpha_T'] else max(1e-15,max(0.,1-p['k_fg']*np.log10(amax/p['alpha_T']))**2)
            row.update(H_max_final=hmax,alpha_max_final=amax,f_min_final=float(fmin))
        table.append(row)
    write_csv(ROOT/"verification/study_summary.csv",table)
    write_csv(ROOT/"verification/output_audit.csv",audits)

    mesh=[];time=[];fields=[];time_fields=[];quad=[];history=[]
    mesh_ref="D0_structured_64_s1_c1_spc80" if a.version=="active_tie" else "D0_structured_64_s1_c5_spc80"
    if mesh_ref in obs:
        for label,rec in records:
            if rec["n"] in (8,16,32,64) and rec["steps_per_cycle"]==80 and label in obs and rec["cycles"] in (1,5):
                limit=80 if a.version=="active_tie" else 400
                cycle=1 if a.version=="active_tie" else 5
                mesh.append(dict(label=label,reference=mesh_ref,**compare_observables(obs[label][:limit],obs[mesh_ref][:limit])))
                if (label,cycle) in states and (mesh_ref,cycle) in states:
                    fields.append(dict(label=label,reference=mesh_ref,cycle=cycle,**compare_fields(states[(label,cycle)],states[(mesh_ref,cycle)])))
    time_ref=next((f'D0_structured_16_s1_c5_spc{s}' for s in (1280,640,320)
                   if f'D0_structured_16_s1_c5_spc{s}' in obs),'D0_structured_16_s1_c5_spc320')
    if time_ref in obs:
        for label,rec in records:
            if rec["cycles"]==5 and rec["mesh"]=="structured" and rec["n"]==16 and label in obs:
                time.append(dict(label=label,reference=time_ref,**compare_observables(obs[label],obs[time_ref])))
                for cyc in (1,5):
                    if (label,cyc) in states and (time_ref,cyc) in states:
                        time_fields.append(dict(label=label,reference=time_ref,cycle=cyc,
                                                **compare_fields(states[(label,cyc)],states[(time_ref,cyc)])))
                        s=states[(label,cyc)]["svars"];r=states[(time_ref,cyc)]["svars"]
                        hh={}
                        state_columns=[('H',hcols),('alpha',acols),('W',hcols+2)]
                        if state_stride==4:
                            state_columns.append(('psi_prev',hcols+3))
                        for key,cols in state_columns:
                            e=s[:,cols]-r[:,cols]
                            hh[key+"_relative_l2"]=float(np.linalg.norm(e)/max(np.linalg.norm(r[:,cols]),1e-30))
                        history.append(dict(label=label,reference=time_ref,cycle=cyc,**hh))
    for n in (16,32):
        label=f"D0_structured_{n}_s1_c5_spc80"
        qpath=ROOT/f"verification/quadrature_comparison/{a.version}/n{n}_3x3/observables.csv"
        if qpath.exists() and label in obs:
            q=np.atleast_1d(np.genfromtxt(qpath,delimiter=",",names=True))
            quad.append(dict(n=n,reference="3x3 quadrature",cycles=3,
                             **compare_observables(obs[label][:len(q)],q)))
    for name,data in [("mesh_convergence",mesh),("time_convergence",time),("field_convergence",fields),
                      ("time_field_convergence",time_fields),
                      ("history_time_convergence",history),("quadrature_comparison",quad)]:
        write_csv(ROOT/f"verification/{name}.csv",data)
    sources=contract['source'] if isinstance(contract['source'],list) else [contract['source']]
    result=dict(version=a.version,
                element_source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                       for name in sources+['element/uel_contract.json']},
                completed_runs=len(records),ok_runs=sum(r["status"]=="ok" for _,r in records),
                output_files_audited=len(audits),mesh_reference=mesh_ref,time_reference=time_ref,
                all_differences_are_self_convergence_not_independent_full_PDE_error=True,
                max_intact_constraint=max((x["intact_constraint_max"] for x in audits),default=None),
                max_bridge_relative_error=max((x["bridge_max_relative_error"] for x in audits),default=None))
    (ROOT/"verification/study_analysis.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    if a.require_complete:
        assert len(records)==20 and result["ok_runs"]==20,"not all planned studies completed"
        assert len(mesh)==8 and len(time)>=8 and len(time_fields)>=16 and len(quad)==2


if __name__=="__main__":
    main()
