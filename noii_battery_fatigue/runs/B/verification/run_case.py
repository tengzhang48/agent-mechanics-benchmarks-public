"""Additional runs on the unmodified host, with detailed acceptance diagnostics.

The supplied host/run_problem.py is used separately for the public runs. This
driver uses the same fehost.solve and ContractUEL and only extends observation.
"""
import argparse
import csv
from dataclasses import asdict
from pathlib import Path
import hashlib
import json
import time

import numpy as np
from scipy.sparse.linalg import spsolve

from common import ROOT, compiled, fehost
import meshes
from observables import observables


def run(problem, mesh="structured", n=8, seed=1, cycles=None, spc=None, out=None,
        uniform_strain=None, tol=1e-8, module_contract=None, source_manifest=None):
    module, contract = compiled() if module_contract is None else module_contract
    if source_manifest is None and module_contract is None:
        names = list(contract['source']) if isinstance(contract['source'], list) else [contract['source']]
        names += ['element/uel_contract.json', 'element/model.py', 'element/generate.py']
        if (ROOT/'element/lump_reaction.py').exists():
            names.append('element/lump_reaction.py')
        source_manifest = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                           for name in names}
    prob = json.loads(Path(problem).read_text()) if not isinstance(problem, dict) else problem
    prm, load = prob["parameters"], prob["loading"]
    contract = {**contract, "props":[prm.get(k,v) for k,v in zip(contract["prop_names"],contract["props"])]}
    cycles = cycles or load["N_cyc"]
    spc = spc or load["steps_per_cycle"]
    dt=load["T"]/spc
    ns=cycles*spc
    nodes, elems = meshes.structured(n) if mesh=="structured" else meshes.distorted(n,seed=seed)
    el=fehost.ContractUEL(module,contract,len(elems),dt)
    dpn=el.dof_per_node
    comp={k:el.comp(k) for k in fehost.REQUIRED_DOFS}
    x,y=nodes.T
    edges={"left":np.isclose(x,0),"right":np.isclose(x,1),"bottom":np.isclose(y,0),"top":np.isclose(y,1)}
    cedges=prob.get("c_edges",["top"])
    csel=np.ones(len(x),bool) if "all_nodes" in cedges else np.zeros(len(x),bool)
    for edge in (edges if "all" in cedges else cedges):
        if edge in edges:
            csel |= edges[edge]
    cnodes=np.where(csel)[0]
    notch=meshes.notch_nodes(nodes,prob["notch"]["x_end"]) if prob.get("notch") else np.array([],int)
    cons=[(meshes.node_at(nodes,xx,yy),cm) for xx,yy,cm in prob.get("constraints",[])]

    def bc(frac):
        t=frac*ns*dt
        cb=prm["c0"]+load["A"]*(1-np.cos(2*np.pi*t/load["T"]))/2
        b={dpn*k+comp["c"]:cb for k in cnodes}
        b.update({dpn*k+comp["d"]:1. for k in notch})
        for k,cm in cons:
            for cname in cm:
                b[dpn*k+comp[cname]]=0.
        if uniform_strain is not None:
            amp=uniform_strain*(1-np.cos(2*np.pi*t/load["T"]))/2
            # Isotropic in-plane tensile loading; all u,c values prescribed.
            for k in range(len(x)):
                b[dpn*k+comp["ux"]]=amp*x[k]
                b[dpn*k+comp["uy"]]=amp*y[k]
                b[dpn*k+comp["c"]]=prm["c0"]
        return b

    U0=np.zeros(len(x)*dpn);U0[comp["c"]::dpn]=prm["c0"]
    U0[dpn*notch+comp["d"]]=1.
    out=Path(out or ROOT/"verification"/"custom_runs"/f'{prob["id"]}_{mesh}_{n}_spc{spc}')
    out.mkdir(parents=True,exist_ok=True)
    rows=[];increments=[];attempt_diagnostics=[]
    counts={"linear_solves":0,"uel_batch_calls":0}
    last_trial={}
    old={"U":U0.copy(),"sv":el.svars.copy()}
    minimum={"H":0.,"alpha":0.,"nodal_damage":0.}
    slots={key:np.array(v,int)-1 for key,v in contract["history_slots"].items()}

    def observe(step,frac,U,metrics=None):
        d=U[comp["d"]::dpn].copy()
        rows.append(dict(step=int(step),t=frac*ns*dt,d_max=float(np.max(d)),d_min=float(np.min(d)),
                         **observables(nodes,elems,U,dpn,comp,prm)))
        if step%spc==0 or step==ns:
            np.savez(out/f'fields_cycle_{-(-step//spc):03d}.npz',nodes=nodes,elems=elems,U=U,
                     dof_names=el.dof_names,svars=el.svars.copy())

    def commit(nodes_,elems_,dpn_,U,previous):
        el.commit_step(nodes_,elems_,dpn_,U,previous)
        for key,s in slots.items():
            minimum[key]=min(minimum[key],float(np.min(el.svars[:,s]-old["sv"][:,s])))
        minimum["nodal_damage"]=min(minimum["nodal_damage"],float(np.min(U[comp["d"]::dpn]-previous[comp["d"]::dpn])))
        increments.append(dict(t=float(el._time[1]),dt=float(el.dt),
                               H_max=float(np.max(el.svars[:,slots["H"]])),
                               alpha_max=float(np.max(el.svars[:,slots["alpha"]])),
                               d_max=float(np.max(U[comp["d"]::dpn])),d_min=float(np.min(U[comp["d"]::dpn])),
                               observed_linear_solves=counts["linear_solves"],
                               observed_uel_batch_calls=counts["uel_batch_calls"],
                               intact_constraint_max=(float(np.max(abs((U[comp["d"]::dpn]-1.)+
                                                          U[el.comp("vfrac")::dpn])))
                                                      if "vfrac" in el.dof_names else None)))
        old["sv"]=el.svars.copy();old["U"]=U.copy()

    def counted_batch(coords_,U_,DU_):
        counts["uel_batch_calls"]+=1
        # Capture input before evaluation, including a PNEWDT rejection. A
        # rejected call otherwise leaves only an earlier successful trial.
        last_trial.update(Ue=U_.copy(),DUe=DU_.copy(),svars=el.svars.copy(),
                          time=el._time.copy(),dt=el.dt,last_call_success=False)
        last_trial.pop('Re',None)
        last_trial.pop('Ke',None)
        result=el.element_rk_batch(coords_,U_,DU_)
        last_trial.update(Re=result[0].copy(),Ke=result[1].copy(),last_call_success=True)
        return result

    def counted_solve(K,b):
        counts["linear_solves"]+=1
        return spsolve(K,b)

    def attempted(rec):
        attempt_diagnostics.append(dict(**asdict(rec),observed_linear_solves=counts["linear_solves"],
                                        observed_uel_batch_calls=counts["uel_batch_calls"]))
        counts.update(linear_solves=0,uel_batch_calls=0)

    status="ok"; report=None;t0=time.monotonic();U=None
    try:
        # Exactly fehost.solve's monolithic setup, with the public runtime's
        # linear_solve/attempt_callback hooks used solely for counting work.
        U,report=fehost.solve_steps_staggered(
            nodes,elems,dpn,el.element_rk,bc,[list(range(dpn))],ns,
            batch_fn=counted_batch,linear_solve=counted_solve,commit_fn=commit,
            state=fehost.CompiledElementState([el]),dt=dt,tol=tol,max_cutbacks=10,
            cutback_on_element_rejection=True,residual_scale=fehost.residual_weights(el,contract),
            U0=U0,step_callback=observe,max_field_newton=25,attempt_callback=attempted,
            return_report=True)
    except Exception as exc:
        status=f"failed: {type(exc).__name__}: {exc}"
        report=getattr(exc,"report",None)
        np.savez(out/"failed_trial.npz",nodes=nodes,elems=elems,**last_trial)
    np.savez(out/"last_state.npz",nodes=nodes,elems=elems,U=old["U"],svars=old["sv"],
             dof_names=el.dof_names)
    summary=dict(status=status,problem=prob,mesh=mesh,n=n,seed=seed,cycles=cycles,steps_per_cycle=spc,
                 dt=dt,accepted_output_steps=len(rows),wall_s=time.monotonic()-t0,
                 history_min_change=minimum,uniform_strain=uniform_strain,tol=tol,
                 element_source_sha256=source_manifest,
                 source_identity_note=('supplied alternate-module identity' if module_contract is not None
                                       else 'production element at compile/load'))
    if report:
        summary.update(report.summary())
        aa=attempt_diagnostics
        (out/"attempts.json").write_text(json.dumps(aa,indent=1)+"\n")
        summary["max_sweeps"]=max((a.sweeps for a in report.attempts),default=0)
        summary["observed_linear_solves_total"]=sum(x["observed_linear_solves"] for x in aa)
        summary["observed_linear_solves_max"]=max((x["observed_linear_solves"] for x in aa),default=0)
        summary["min_dt"]=min((a.dt for a in report.attempts),default=None)
    for name,data in [("observables",rows),("increments",increments)]:
        if data:
            with (out/(name+".csv")).open("w") as f:
                writer=csv.DictWriter(f,fieldnames=list(data[0]));writer.writeheader();writer.writerows(data)
    (out/"run.json").write_text(json.dumps(summary,indent=1)+"\n")
    print(json.dumps({k:v for k,v in summary.items() if k!="problem"}),flush=True)
    return summary


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--problem",default=str(ROOT/"host/dev_problem/d0.json"))
    ap.add_argument("--mesh",default="structured",choices=["structured","distorted"])
    ap.add_argument("--n",type=int,default=8)
    ap.add_argument("--seed",type=int,default=1)
    ap.add_argument("--cycles",type=int)
    ap.add_argument("--spc",type=int)
    ap.add_argument("--out")
    ap.add_argument("--uniform-strain",type=float)
    ap.add_argument("--tol",type=float,default=1e-8)
    a=ap.parse_args()
    run(a.problem,a.mesh,a.n,a.seed,a.cycles,a.spc,a.out,a.uniform_strain,a.tol)
