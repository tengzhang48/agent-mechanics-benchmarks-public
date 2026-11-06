"""Run a v2 problem with a workspace's UEL on the host.

usage: python run_problem.py WORKSPACE --problem dev_problem/d0.json --mesh structured|distorted
                             --n 16 [--seed 0] [--cycles K] [--steps-per-cycle 40] [--out DIR]

The problem file gives the parameters, the loading c_b(t) = c0 + A[1 - cos(2 pi t / T)]/2, the edges
that carry c_b ("left", "right", "bottom", "top" or "all"), an optional notch (d = 1 on a segment of
y = const), and the displacement constraints [[x, y, ["ux", "uy"]], ...] at mesh nodes.
Writes DIR/observables.csv (one row per accepted output step), DIR/fields_cycle_XXX.npz and DIR/run.json.
Observables are computed by the host from nodal values with Q4 interpolation.
"""
import argparse
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fehost  # noqa: E402
import meshes  # noqa: E402
from observables import observables  # noqa: E402


def edge_nodes(nodes, edges):
    x, y = nodes[:, 0], nodes[:, 1]
    test = {"left": np.isclose(x, 0), "right": np.isclose(x, 1),
            "bottom": np.isclose(y, 0), "top": np.isclose(y, 1)}
    if "all" in edges:
        edges = list(test)
    mask = np.zeros(len(nodes), bool)
    for e in edges:
        mask |= test[e]
    return np.where(mask)[0]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("workspace")
    ap.add_argument("--problem", default=os.path.join(HERE, "dev_problem", "d0.json"))
    ap.add_argument("--mesh", default="structured", choices=["structured", "distorted"])
    ap.add_argument("--n", type=int, default=16)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--cycles", type=int, default=None)
    ap.add_argument("--steps-per-cycle", type=int, default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)

    prob = json.load(open(a.problem))
    prm, load = prob["parameters"], prob["loading"]
    ncyc = a.cycles or load["N_cyc"]
    spc = a.steps_per_cycle or load["steps_per_cycle"]
    dt = load["T"] / spc
    name = prob.get("id", "problem")
    out = a.out or os.path.join(a.workspace, "verification", "runs",
                                f"{name}_{a.mesh}_{a.n}_s{a.seed}_spc{spc}")
    os.makedirs(out, exist_ok=True)

    contract = fehost.load_contract(a.workspace)
    props = list(contract["props"])
    for i, pname in enumerate(contract["prop_names"]):
        if pname in prm:
            props[i] = prm[pname]
    contract = {**contract, "props": props}
    nodes, elems = (meshes.structured(a.n) if a.mesh == "structured"
                    else meshes.distorted(a.n, seed=a.seed))
    mod = fehost.build_module(a.workspace, contract, os.path.join(out, "build"),
                              f"uel_{name}_{a.mesh}_{a.n}_{a.seed}_{os.getpid()}")
    el = fehost.ContractUEL(mod, contract, len(elems), dt)
    dpn = el.dof_per_node
    comp = {k: el.comp(k) for k in fehost.REQUIRED_DOFS}
    cnodes = edge_nodes(nodes, prob["c_edges"])
    notch = (meshes.notch_nodes(nodes, prob["notch"]["x_end"]) if prob.get("notch")
             else np.array([], dtype=int))
    cons = [(meshes.node_at(nodes, x, y), comps) for x, y, comps in prob["constraints"]]
    c0, A, T = prm["c0"], load["A"], load["T"]
    n_steps = ncyc * spc

    def dirichlet(frac):
        t = frac * n_steps * dt
        cb = c0 + A * (1 - np.cos(2 * np.pi * t / T)) / 2
        bc = {dpn * k + comp["c"]: cb for k in cnodes}
        bc.update({dpn * k + comp["d"]: 1.0 for k in notch})
        for k, comps in cons:
            for cname in comps:
                bc[dpn * k + comp[cname]] = 0.0
        return bc

    U0 = np.zeros(len(nodes) * dpn)
    U0[comp["c"]::dpn] = c0
    U0[dpn * notch + comp["d"]] = 1.0
    rows = []
    t0 = time.time()
    slots = {k: np.asarray(v, int) - 1 for k, v in (contract.get("history_slots") or {}).items()}
    hist = {"prev": {k: el.svars[:, s].copy() for k, s in slots.items()},
            "min_change": {k: np.inf for k in slots}}

    def step_cb(step, frac, U, metrics=None):
        d = U[comp["d"]::dpn].copy()
        if len(notch):
            d[notch] = -np.inf
        for k, s in slots.items():   # history must never decrease between accepted steps
            now = el.svars[:, s]
            hist["min_change"][k] = min(hist["min_change"][k], float((now - hist["prev"][k]).min()))
            hist["prev"][k] = now.copy()
        rows.append({"step": int(step), "t": frac * n_steps * dt, "d_max": float(d.max()),
                     **observables(nodes, elems, U, dpn, comp, prm)})
        if step % spc == 0 or step == n_steps:  # the driver numbers steps from 1
            np.savez(os.path.join(out, f"fields_cycle_{-(-step // spc):03d}.npz"),
                     nodes=nodes, elems=elems, U=U, dof_names=np.array(el.dof_names),
                     svars=el.svars.copy())

    status, report = "ok", None
    try:
        U, report = fehost.solve(nodes, elems, el, dirichlet, n_steps, dt, U0=U0,
                                 residual_scale=fehost.residual_weights(el, contract),
                                 step_callback=step_cb)
    except Exception as exc:  # record the failure; the observables so far are kept
        status = f"failed: {type(exc).__name__}: {exc}"[:2000]
    keys = ["step", "t", "d_max", "A_damage", "D_int", "C_int", "E_el"]
    with open(os.path.join(out, "observables.csv"), "w") as f:
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(f"{r[k]:.12g}" for k in keys) + "\n")
    summary = {"status": status, "problem": name, "mesh": a.mesh, "n": a.n, "seed": a.seed,
               "cycles": ncyc, "steps_per_cycle": spc, "dt": dt,
               "accepted_output_steps": len(rows), "wall_s": round(time.time() - t0, 2),
               "history_min_change": {k: (None if not np.isfinite(v) else v)
                                      for k, v in hist["min_change"].items()}}
    if report is not None:
        summary.update(report.summary())
        summary["rejected_attempts"] = sum(1 for x in report.attempts if not x.converged)
        summary["min_dt"] = min((x.dt for x in report.attempts), default=None)
    json.dump(summary, open(os.path.join(out, "run.json"), "w"), indent=1)
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
