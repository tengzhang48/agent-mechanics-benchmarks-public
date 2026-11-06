"""Global studies on the unmodified host, with a reusable compiled UEL.

The solve is fehost.solve (one monolithic field group), with its default
tolerances, line search and cutbacks. Only reporting is added here.
"""
from pathlib import Path
import sys
import json
import argparse
import time
import hashlib
from dataclasses import asdict
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from verification.common import compiled
import fehost
import meshes
from observables import observables
from run_problem import edge_nodes


def run_case(n=16, mesh="structured", seed=0, cycles=5, spc=80,
             problem=None, out=None, solver_stats=False, contractfile=None):
    problem = Path(problem) if problem else ROOT/"host/dev_problem/d0.json"
    prob = json.loads(problem.read_text())
    prm, load = prob["parameters"], prob["loading"]
    c = json.loads(Path(contractfile).read_text()) if contractfile else fehost.load_contract(str(ROOT))
    tag=c.get("formulation_id","")
    prefix=f"{prob['id']}_{tag}_" if tag else f"{prob['id']}_"
    out = Path(out) if out else ROOT/"verification/runs"/f"{prefix}{mesh}_{n}_s{seed}_spc{spc}_cy{cycles}"
    out.mkdir(parents=True,exist_ok=True)
    c = {**c, "props": [prm.get(k, v) for k,v in zip(c["prop_names"],c["props"])]}
    originals=list(c["source"])
    snapshot=out/"source_snapshot";snapshot.mkdir(exist_ok=True)
    snapshots=[]
    for i,s in enumerate(originals):
        target=snapshot/(str(i)+"_"+Path(s).name)
        target.write_bytes((ROOT/s).read_bytes())
        snapshots.append(str(target.resolve().relative_to(ROOT)))
    c={**c,"source":snapshots}
    source_hashes={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in c["source"]}
    module = compiled(c)
    nodes, elems = meshes.structured(n) if mesh=="structured" else meshes.distorted(n,seed)
    dt = load["T"]/spc
    el = fehost.ContractUEL(module,c,len(elems),dt)
    dpn = el.dof_per_node
    comp = {k:el.comp(k) for k in fehost.REQUIRED_DOFS}
    cnodes = edge_nodes(nodes,prob["c_edges"])
    notch = meshes.notch_nodes(nodes,prob["notch"]["x_end"]) if prob.get("notch") else np.array([],int)
    constraints = [(meshes.node_at(nodes,x,y),comps) for x,y,comps in prob["constraints"]]
    def bc(frac):
        t = frac*cycles*load["T"]
        cb = prm["c0"]+load["A"]*(1-np.cos(2*np.pi*t/load["T"]))/2
        result = {dpn*k+comp["c"]:cb for k in cnodes}
        result.update({dpn*k+comp["d"]:1. for k in notch})
        for k, comps in constraints:
            result.update({dpn*k+comp[name]:0. for name in comps})
        return result
    U0 = np.zeros(len(nodes)*dpn)
    U0[comp["c"]::dpn] = prm["c0"]
    U0[dpn*notch+comp["d"]] = 1.
    histories = {k:np.array(v,int)-1 for k,v in c["history_slots"].items()}
    minima = {k:np.inf for k in histories}
    rows = []
    commit_count = 0
    max_d = 0.
    min_d = np.inf
    def commit(nodes_,elems_,dpn_,U,Uold):
        nonlocal commit_count
        old = el.svars.copy()
        el.commit_step(nodes_,elems_,dpn_,U,Uold)
        commit_count += 1
        for key, idx in histories.items():
            delta = el.svars[:,idx]-old[:,idx]
            minima[key] = min(minima[key],float(delta.min()))
            assert np.all(delta>=0), (key,delta.min())
    def callback(step,frac,U,metrics=None):
        nonlocal max_d,min_d
        d = U[comp["d"]::dpn].copy()
        if len(notch):
            d[notch] = -np.inf
        max_d = max(max_d,float(d.max()))
        min_d = min(min_d,float(U[comp["d"]::dpn].min()))
        row = dict(step=step,t=frac*cycles*load["T"],d_max=float(d.max()),
                   **observables(nodes,elems,U,dpn,comp,prm))
        rows.append(row)
        if step%spc==0:
            np.savez_compressed(out/f"fields_cycle_{step//spc:03d}.npz",nodes=nodes,elems=elems,
                                U=U,dof_names=np.array(el.dof_names),svars=el.svars)
        if step%spc==0:
            print(f"{out.name}: cycle {step//spc}, dmax={row['d_max']:.7g}, "
                  f"D={row['D_int']:.7g}, C={row['C_int']:.7g}",flush=True)
    start = time.monotonic()
    status = "ok"; report = None
    newton_stats=[]
    last_trial={}
    try:
        if solver_stats:
            # Exactly fehost.solve's arguments, plus passive counters around
            # the same linear solve and batch calls. No host modification.
            from hostfe.driver import _scipy_solve
            count={"linear_solves":0,"assemblies":0}
            def linear(K,b):
                count["linear_solves"]+=1
                return _scipy_solve(K,b)
            def batch(coords,U,DU):
                count["assemblies"]+=1
                last_trial.update(coords=coords.copy(),U=U.copy(),DU=DU.copy(),
                    svars=el.svars.copy(),dt=np.array(el.dt),time=el._time.copy(),props=el.props.copy())
                return el.element_rk_batch(coords,U,DU)
            def attempt(rec):
                newton_stats.append(dict(step=rec.step,dt=rec.dt,
                    converged=rec.converged,**count))
                count.update(linear_solves=0,assemblies=0)
            U,report=fehost.solve_steps_staggered(nodes,elems,dpn,el.element_rk,bc,
                [list(range(dpn))],cycles*spc,batch_fn=batch,commit_fn=commit,
                state=fehost.CompiledElementState([el]),dt=dt,tol=1e-8,max_cutbacks=10,
                cutback_on_element_rejection=True,residual_scale=fehost.residual_weights(el,c),
                U0=U0,step_callback=callback,max_field_newton=25,return_report=True,
                linear_solve=linear,attempt_callback=attempt)
        else:
            U,report = fehost.solve(nodes,elems,el,bc,cycles*spc,dt,U0=U0,
                residual_scale=fehost.residual_weights(el,c),step_callback=callback,commit_fn=commit)
    except Exception as exc:
        status = f"failed: {type(exc).__name__}: {exc}"
        report = getattr(exc,"report",None)
        if last_trial:
            np.savez_compressed(out/"failed_trial.npz",**last_trial)
    keys = ["step","t","d_max","A_damage","D_int","C_int","E_el"]
    with (out/"observables.csv").open("w") as f:
        f.write(",".join(keys)+"\n")
        for row in rows:
            f.write(",".join(f"{row[k]:.16g}" for k in keys)+"\n")
    summary = dict(status=status,problem=prob["id"],mesh=mesh,n=n,seed=seed,cycles=cycles,
        steps_per_cycle=spc,dt=dt,wall_s=time.monotonic()-start,accepted_output_steps=len(rows),
        committed_subincrements=commit_count,history_min_change=minima,
        d_max_over_run=max_d,d_min_over_run=min_d,
        parameters=prm,host="unmodified fehost.solve, default tol=1e-8",
        source_sha256=source_hashes,compiled_module=module.__name__,original_source_paths=originals)
    (out/"contract.json").write_text(json.dumps(c,indent=2)+"\n")
    if newton_stats:
        (out/"newton_counts.json").write_text(json.dumps(newton_stats,indent=1)+"\n")
        accepted_counts=[r["linear_solves"] for r in newton_stats if r["converged"]]
        summary["newton_linear_solves"]={"min":min(accepted_counts),"max":max(accepted_counts),
            "mean":float(np.mean(accepted_counts)),"total":sum(accepted_counts)}
    if report:
        summary.update(report.summary())
        attempts = report.attempts
        summary.update(rejected_attempts=sum(not a.converged for a in attempts),
                       min_dt=min(a.dt for a in attempts),
                       max_sweeps=max(a.sweeps for a in attempts),
                       max_accepted_sweeps=max(a.sweeps for a in attempts if a.converged))
        with (out/"solver_report.json").open("w") as f:
            json.dump(asdict(report),f,indent=1)
    (out/"run.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary),flush=True)
    return summary


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--n",type=int,default=16)
    ap.add_argument("--mesh",choices=["structured","distorted"],default="structured")
    ap.add_argument("--seed",type=int,default=0)
    ap.add_argument("--cycles",type=int,default=5)
    ap.add_argument("--spc",type=int,default=80)
    ap.add_argument("--problem")
    ap.add_argument("--out")
    ap.add_argument("--solver-stats",action="store_true")
    ap.add_argument("--contract",dest="contractfile")
    args=ap.parse_args()
    r=run_case(**vars(args))
    if r["status"]!="ok":
        raise SystemExit(1)
