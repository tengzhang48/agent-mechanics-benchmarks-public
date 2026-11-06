"""Robustness battery of the battery benchmark (test 1), on the public abaqus_ufl.fe host.

usage (each test is one process, so a Fortran ERROR STOP cannot stop the scorer):
  python battery.py hidden_state WORKSPACE OUT.json
  python battery.py dtime_zero   WORKSPACE OUT.json

The problem is the development problem D0 (smooth, no crack), or the file in V2_PROBLEM (score.py sets
it). The official scorer's battery.py with one change: the solver and the rejection signal are imported
from the public abaqus_ufl.fe instead of the host's vendored copy.

hidden_state: the problem on an 8x8 mesh, 3 increments of the benchmark step. Two runs both reject increment
  2 and retry it at half the step: "after" rejects in commit, after the UEL saw that attempt's trial
  states; "before" rejects at the attempt's first assembly, before the UEL is called. A UEL whose
  only memory is SVARS gives identical U and SVARS.
dtime_zero: one direct call with DTIME = 0 at the state after 3 increments; it must return
  finite RHS/AMATRX or PNEWDT < 1, and must not crash.
"""
import json
import os
import sys

import numpy as np

HOST = os.environ.get("V2_HOST", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "host"))
sys.path.insert(0, HOST)
import fehost  # noqa: E402
import meshes  # noqa: E402
from abaqus_ufl.fe import IncrementTooLarge, solve_steps_staggered  # noqa: E402
from run_problem import edge_nodes  # noqa: E402

DEV = json.load(open(os.environ.get("V2_PROBLEM", os.path.join(HOST, "dev_problem", "d0.json"))))


def setup(workspace, n=8, tag="b"):
    contract = fehost.load_contract(workspace)
    prm = DEV["parameters"]
    props = list(contract["props"])
    for i, name in enumerate(contract["prop_names"]):
        if name in prm:
            props[i] = prm[name]
    contract = {**contract, "props": props}
    nodes, elems = meshes.structured(n)
    work = os.path.join(os.path.dirname(os.path.abspath(sys.argv[3])), f"build_{tag}")
    mod = fehost.build_module(workspace, contract, work, f"bat_{tag}_{os.getpid()}")
    return contract, nodes, elems, mod


def d1_bc(nodes, el, t):
    dpn = el.dof_per_node
    comp = {k: el.comp(k) for k in fehost.REQUIRED_DOFS}
    c0, A, T = DEV["parameters"]["c0"], DEV["loading"]["A"], DEV["loading"]["T"]
    cb = c0 + A * (1 - np.cos(2 * np.pi * t / T)) / 2
    bc = {dpn * k + comp["c"]: cb for k in edge_nodes(nodes, DEV["c_edges"])}
    bc.update({dpn * k + comp["d"]: 1.0 for k in notch_of(nodes)})
    for x, y, comps in DEV["constraints"]:
        k = meshes.node_at(nodes, x, y)
        for cname in comps:
            bc[dpn * k + comp[cname]] = 0.0
    return bc


def notch_of(nodes):
    if DEV.get("notch"):
        return meshes.notch_nodes(nodes, DEV["notch"]["x_end"])
    return np.array([], dtype=int)


def initial(nodes, el):
    dpn = el.dof_per_node
    U0 = np.zeros(len(nodes) * dpn)
    U0[el.comp("c")::dpn] = DEV["parameters"]["c0"]
    U0[dpn * notch_of(nodes) + el.comp("d")] = 1.0
    return U0


def run_three(contract, nodes, elems, mod, mode):
    dt = DEV["loading"]["T"] / DEV["loading"]["steps_per_cycle"]
    el = fehost.ContractUEL(mod, contract, len(elems), dt)
    st = {"commits": 0, "armed": False, "fired": False}

    def commit_fn(nodes_, elems_, dpn, U, U_prev):
        st["commits"] += 1
        if mode == "after" and st["commits"] == 2 and not st["fired"]:
            st["fired"] = True
            raise IncrementTooLarge("forced rejection after the UEL saw the attempt")
        el.commit_step(nodes_, elems_, dpn, U, U_prev)
        if mode == "before" and st["commits"] == 1:
            st["armed"] = True

    def batch_fn(coords_all, U_all, DU_all):
        if st["armed"] and not st["fired"]:
            st["fired"] = True
            raise IncrementTooLarge("forced rejection before the UEL saw the attempt")
        return el.element_rk_batch(coords_all, U_all, DU_all)

    n_steps = 3
    U, rep = solve_steps_staggered(
        nodes, elems, el.dof_per_node, el.element_rk, lambda f: d1_bc(nodes, el, f * n_steps * dt),
        [list(range(el.dof_per_node))], n_steps, batch_fn=batch_fn, commit_fn=commit_fn,
        state=fehost.CompiledElementState([el]), dt=dt, tol=1e-8, max_cutbacks=4,
        cutback_on_element_rejection=True, U0=initial(nodes, el),
        residual_scale=fehost.residual_weights(el, contract), return_report=True)
    return U, el.svars.copy(), st["fired"], el


def hidden_state(workspace, out):
    contract, nodes, elems, mod = setup(workspace, tag="hs")
    try:
        Ua, sva, fa, _ = run_three(contract, nodes, elems, mod, "after")
        Ub, svb, fb, _ = run_three(contract, nodes, elems, mod, "before")
    except Exception as exc:  # a solve that fails here is a failed test, with its reason
        res = {"test": "hidden_state", "pass": False,
               "error": f"{type(exc).__name__}: {exc}"[:1000]}
        json.dump(res, open(out, "w"), indent=1)
        print(json.dumps(res))
        return
    res = {"test": "hidden_state", "fired": [fa, fb],
           "max_abs_dU": float(np.abs(Ua - Ub).max()),
           "max_abs_dSVARS": float(np.abs(sva - svb).max())}
    res["pass"] = bool(fa and fb and res["max_abs_dU"] <= 1e-12 and res["max_abs_dSVARS"] <= 1e-12)
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res))


def dtime_zero(workspace, out):
    contract, nodes, elems, mod = setup(workspace, tag="dz")
    U, sv, _, el = run_three(contract, nodes, elems, mod, None)
    dpn = el.dof_per_node
    e = int(np.argmax([np.abs(U[el.comp("d")::dpn][en]).sum() for en in elems]))  # a damaged element
    gm = (elems[e][:, None] * dpn + np.arange(dpn)[None, :]).ravel()
    el.dt = 0.0
    rhs, amatrx, _, pn = el._call(e, nodes[elems[e]], U[gm], np.zeros(len(gm)), sv[e].copy())
    finite = bool(np.all(np.isfinite(rhs)) and np.all(np.isfinite(amatrx)))
    res = {"test": "dtime_zero", "element": e, "finite": finite, "pnewdt": float(pn),
           "pass": bool(finite or pn < 1.0)}
    json.dump(res, open(out, "w"), indent=1)
    print(json.dumps(res))


if __name__ == "__main__":
    {"hidden_state": hidden_state, "dtime_zero": dtime_zero}[sys.argv[1]](sys.argv[2], sys.argv[3])
