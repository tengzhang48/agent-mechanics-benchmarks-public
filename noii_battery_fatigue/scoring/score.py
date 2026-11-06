"""The element scorer of the battery benchmark (test 1), for any problem file and reference folder.

usage:
  python score.py run   WORKSPACE OUT --problem PROBLEM.json --key KEYDIR [--stress-problem P.json] [--jobs 8]
  python score.py score WORKSPACE OUT --problem PROBLEM.json --key KEYDIR [--stress-problem P.json]

This is the scorer that produced the official scores, with three changes and no change to any check or
tolerance: the problem file, the reference folder and the optional stress variant are arguments (they were
fixed paths to the hidden problem); the host is ../host (the public abaqus_ufl.fe host); the Python is the
one running this script. The released practice run is the development problem D0 with its reference
(--problem ../host/dev_problem/d0.json --key reference_D0). The stress variant belongs to the hidden problem
and is not released, so on D0 the check R5 (damage near 1) is reported as not run.

Runs (PROBLEM, its first 10 cycles):
  s16_40, s16_80, s32_40, s32_80   structured n = 16, 32 at T/40 and T/80
  dist0, dist1                     distorted n = 32 (seeds 0, 1), T/40
  big                              structured n = 32, T/10 (a quarter of the steps)
  stress                           only with --stress-problem: n = 32, T/40, 2 cycles; d must reach 0.9
  battery                          hidden_state and dtime_zero on PROBLEM (8x8)

The key is the dolfinx reference F1 on u64 (KEYDIR/F1_u64_div{40,80}); its spatial error is ~1e-4.
The element's history update is its own choice, so at a given dt it need not match the key's lagged
scheme. Both are therefore compared in the dt limit, v_ext = 2 v(T/80) - v(T/40), at the cycle ends.
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
HOST = os.environ.get("V2_HOST", os.path.join(os.path.dirname(HERE), "host"))
PY = sys.executable
PROB = PROB_STRESS = KEY = None   # set from the command line in main()
FIELDS = ("ux", "uy", "c", "d")
OBS = ("d_max", "D_int", "C_int", "E_el")
NCYC = 10

# Tolerances, set before the official runs from the key's errors and a pipeline check with an earlier element:
# a consistent Q4 element (v1 at chi = 0) is within 1.3e-5 in C_int and 1e-4 RMS in u, c at n = 32, but 8 %
# off in E_el (a small remainder of nearly cancelling strains), so E_el is reported, not scored.
# v1 at chi = 0.25 (incomplete grad p_eff) is 1.5 % off in C_int and 1.7-3.8 % RMS in c, u.
TOL = {
    "obs": {"d_max": 0.04, "D_int": 0.04, "C_int": 0.005},  # max over cycle ends / max |key|
    "field_rms": 0.01,          # RMS of the nodal difference / max |key|, max over cycle ends
    "drop": 1.5,                # C2: the error against the key must drop at least this much from n = 16 to 32 ...
    "converged": 0.005,         # ... unless it is already below this at n = 32 (above the key's own residual)
    "history": -1e-12,          # R1: smallest allowed change of a history slot between accepted steps
    "d_range": (-0.01, 1.01),   # R5: d stays in this range ...
    "d_reached": 0.9,           # ... and reaches at least this
}

RUNS = {
    "s16_40": dict(n=16, spc=40), "s16_80": dict(n=16, spc=80),
    "s32_40": dict(n=32, spc=40), "s32_80": dict(n=32, spc=80),
    "dist0": dict(n=32, spc=40, mesh="distorted", seed=0),
    "dist1": dict(n=32, spc=40, mesh="distorted", seed=1),
    "big": dict(n=32, spc=10),
    "stress": dict(n=32, spc=40, problem=None, cycles=2),   # problem: --stress-problem
}


# ------------------------------------------------------------------------------------------------ run
def commands(ws, out):
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "PYTHONDONTWRITEBYTECODE": "1",
           "V2_PROBLEM": PROB, "V2_HOST": HOST}
    cmds = {}
    for name, r in RUNS.items():
        if name == "stress":
            if PROB_STRESS is None:
                continue
            r = {**r, "problem": PROB_STRESS}
        cmds[name] = [PY, os.path.join(HOST, "run_problem.py"), ws, "--problem", r.get("problem") or PROB,
                      "--mesh", r.get("mesh", "structured"), "--n", str(r["n"]), "--seed", str(r.get("seed", 0)),
                      "--steps-per-cycle", str(r["spc"]), "--cycles", str(r.get("cycles", NCYC)),
                      "--out", os.path.join(out, name)]
    for test in ("hidden_state", "dtime_zero"):
        cmds[f"battery_{test}"] = [PY, os.path.join(HERE, "battery.py"), test, ws,
                                   os.path.join(out, "battery", f"{test}.json")]
    return cmds, env


def run(ws, out, jobs, only=None, timeout_h=2.0):
    os.makedirs(os.path.join(out, "logs"), exist_ok=True)
    os.makedirs(os.path.join(out, "battery"), exist_ok=True)
    cmds, env = commands(ws, out)
    todo = [k for k in cmds if not only or k in only]
    active, done = {}, {}
    while todo or active:
        while todo and len(active) < jobs:
            k = todo.pop(0)
            log = open(os.path.join(out, "logs", f"{k}.log"), "w")
            active[k] = (subprocess.Popen(cmds[k], env=env, stdout=log, stderr=subprocess.STDOUT), time.time(), log)
            print(f"{time.strftime('%H:%M:%S')} start {k}", flush=True)
        for k, (p, t0, log) in list(active.items()):
            rc = p.poll()
            if rc is None and time.time() - t0 > timeout_h * 3600:
                p.kill()
                rc = "timeout"
            if rc is not None:
                log.close()
                done[k] = {"returncode": rc, "wall_s": round(time.time() - t0, 1)}
                del active[k]
                print(f"{time.strftime('%H:%M:%S')} end {k}: {done[k]}", flush=True)
        time.sleep(2)
    json.dump({"workspace": ws, "runs": done, "commands": cmds}, open(os.path.join(out, "runs.json"), "w"), indent=1)


# ------------------------------------------------------------------------------------------------ load
def read_csv(path):
    with open(path) as f:
        keys = f.readline().strip().split(",")
        rows = [line.strip().split(",") for line in f if line.strip()]
    return {k: np.array([float(r[i]) for r in rows]) for i, k in enumerate(keys)}


def cycle_index(t, T, ncyc):
    """Row index of each cycle end k T, k = 1..ncyc (-1 if missing)."""
    return np.array([int(np.argmin(np.abs(t - k * T))) if len(t) and np.min(np.abs(t - k * T)) < 1e-6 * T
                     else -1 for k in range(1, ncyc + 1)])


def load_run(out, name, T, ncyc):
    d = os.path.join(out, name)
    if not os.path.exists(os.path.join(d, "run.json")):
        return None
    r = {"run": json.load(open(os.path.join(d, "run.json"))), "obs": None, "fields": None}
    if os.path.exists(os.path.join(d, "observables.csv")):
        o = read_csv(os.path.join(d, "observables.csv"))
        if "t" in o and len(o["t"]):
            idx = cycle_index(o["t"], T, ncyc)
            r["obs"] = {k: np.where(idx >= 0, o[k][np.maximum(idx, 0)], np.nan) for k in OBS}
    fs = []
    for k in range(1, ncyc + 1):
        p = os.path.join(d, f"fields_cycle_{k:03d}.npz")
        if not os.path.exists(p):
            break
        z = np.load(p)
        names = [str(s) for s in z["dof_names"]]
        U = z["U"].reshape(-1, len(names))
        fs.append(np.stack([U[:, names.index(f)] for f in FIELDS]))
    if fs:
        r["fields"] = np.array(fs)  # (cycles done, 4, N)
    return r


def load_key(div, n_host, T, ncyc):
    d = os.path.join(KEY, f"F1_u64_div{div}")
    o = read_csv(os.path.join(d, "observables.csv"))
    idx = cycle_index(o["t"], T, ncyc)
    obs = {k: np.where(idx >= 0, o[k][np.maximum(idx, 0)], np.nan) for k in ("D_int", "C_int", "E_el")}
    obs["d_max"] = np.where(idx >= 0, o[f"d_max_host{n_host}"][np.maximum(idx, 0)], np.nan)
    z = np.load(os.path.join(d, "host_nodes.npz"))
    fidx = cycle_index(z["t"], T, ncyc)
    fields = z[f"n{n_host}_values"][np.maximum(fidx, 0)].astype(float)  # (ncyc, 4, N): ux, uy, c, d
    fields[fidx < 0] = np.nan
    return {"obs": obs, "fields": fields}


# ------------------------------------------------------------------------------------------------ metrics
def rel_series(a, r):
    """max |a - r| over the cycle ends / max |r|; inf if any compared value is missing or not finite."""
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(r))):
        return float("inf")
    diff, scale = np.max(np.abs(a - r)), np.max(np.abs(r))
    return float(diff / scale) if scale > 0 else (0.0 if diff == 0 else float("inf"))


def field_rms(a, r):
    """Per field: max over cycle ends of RMS(a - r), divided by max |r|."""
    if a is None or a.shape != r.shape or not np.all(np.isfinite(r)):
        return {f: float("inf") for f in FIELDS}
    out = {}
    for i, f in enumerate(FIELDS):
        scale = np.max(np.abs(r[:, i]))
        out[f] = float(np.max(np.sqrt(np.mean((a[:, i] - r[:, i]) ** 2, axis=1))) / scale) if scale > 0 else 0.0
    return out


def complete(r, ncyc, spc):
    return (r is not None and r["run"].get("status") == "ok"
            and r["run"].get("accepted_output_steps") == ncyc * spc and r["obs"] is not None
            and all(np.all(np.isfinite(v)) for v in r["obs"].values()))


def vs_key(R, ok, n, T):
    """Errors of the element's dt limit against the key's dt limit, on mesh n (None if a run is missing)."""
    a40, a80 = R.get(f"s{n}_40"), R.get(f"s{n}_80")
    if not (ok[f"s{n}_40"] and ok[f"s{n}_80"]):
        return None
    k40, k80 = load_key(40, n, T, NCYC), load_key(80, n, T, NCYC)
    obs = {q: rel_series(2 * a80["obs"][q] - a40["obs"][q], 2 * k80["obs"][q] - k40["obs"][q]) for q in OBS}
    fa = None if a40["fields"] is None or a80["fields"] is None else 2 * a80["fields"] - a40["fields"]
    return {"obs": obs, "fields": field_rms(fa, 2 * k80["fields"] - k40["fields"]),
            "same_dt_T40": {q: rel_series(a40["obs"][q], k40["obs"][q]) for q in OBS}}  # info: the key's scheme?


# ------------------------------------------------------------------------------------------------ score
def score(ws, out):
    T = json.load(open(PROB))["loading"]["T"]
    R = {k: load_run(out, k, T, NCYC) for k in RUNS if k != "stress"}
    ok = {k: complete(R[k], NCYC, RUNS[k]["spc"]) for k in R}
    res = {"workspace": os.path.basename(ws), "problem": json.load(open(PROB)).get("id", "problem"), "cycles": NCYC, "tolerances": TOL, "checks": {}, "info": {}}
    C = res["checks"]
    C["R0_runs_complete"] = {"runs": {k: (R[k]["run"]["status"] if R[k] else "missing") for k in R},
                             "pass": all(ok.values())}

    # C1: the dt limit on n = 32 against the key's dt limit
    e = {n: vs_key(R, ok, n, T) for n in (16, 32)}
    if e[32]:
        C["C1_key_observables"] = {"n32": {q: e[32]["obs"][q] for q in OBS},
                                   "pass": all(e[32]["obs"][q] <= TOL["obs"][q] for q in TOL["obs"])}
        C["C1_key_fields"] = {"n32": e[32]["fields"],
                              "pass": all(v <= TOL["field_rms"] for v in e[32]["fields"].values())}
        res["info"]["E_el_vs_key_n32"] = e[32]["obs"]["E_el"]
        res["info"]["same_dt_T40_vs_key_n32"] = e[32]["same_dt_T40"]
    else:
        C["C1_key_observables"] = C["C1_key_fields"] = {"pass": False, "reason": "s32_40 or s32_80 incomplete"}

    # C2: the error against the key drops from n = 16 to n = 32
    if e[16] and e[32]:
        c2 = {}
        for q in TOL["obs"]:
            c2[q] = (e[16]["obs"][q], e[32]["obs"][q])
        for f in FIELDS:
            c2[f"field_{f}"] = (e[16]["fields"][f], e[32]["fields"][f])
        c2 = {k: {"n16": a, "n32": b, "pass": bool(b <= TOL["converged"] or (b > 0 and a / b >= TOL["drop"]))}
              for k, (a, b) in c2.items()}
        C["C2_converges_to_key"] = {"per_quantity": c2, "pass": all(v["pass"] for v in c2.values())}
    else:
        C["C2_converges_to_key"] = {"pass": False, "reason": "a structured run is incomplete"}

    # C3: distortion (distorted n = 32 vs structured n = 32, same dt)
    c3 = {}
    for nm in ("dist0", "dist1"):
        if ok[nm] and ok["s32_40"]:
            c3[nm] = {q: rel_series(R[nm]["obs"][q], R["s32_40"]["obs"][q]) for q in OBS}
            c3[nm]["pass"] = all(c3[nm][q] <= TOL["obs"][q] for q in TOL["obs"])  # E_el reported only
        else:
            c3[nm] = {"pass": False, "reason": "incomplete"}
    C["C3_distortion"] = {"per_mesh": c3, "pass": all(v["pass"] for v in c3.values())}

    # R1: history never decreases (declared slots H and alpha, every run)
    contract = json.load(open(os.path.join(ws, "element", "uel_contract.json")))
    slots = contract.get("history_slots") or {}
    declared = all(k in slots and len(slots[k]) > 0 for k in ("H", "alpha"))
    rs = load_run(out, "stress", T, RUNS["stress"]["cycles"]) if PROB_STRESS else None
    allr = {**R, "stress": rs} if PROB_STRESS else dict(R)
    mins = {k: (v["run"].get("history_min_change") if v else None) for k, v in allr.items()}
    worst = min((x for m in mins.values() if m for x in m.values() if x is not None), default=None)
    C["R1_history"] = {"declared": declared, "slots": slots, "min_change_per_run": mins, "worst": worst,
                       "pass": bool(declared and worst is not None and worst >= TOL["history"])}

    # R2, R3: battery
    for key, test in (("R2_hidden_state", "hidden_state"), ("R3_dtime_zero", "dtime_zero")):
        p = os.path.join(out, "battery", f"{test}.json")
        C[key] = json.load(open(p)) if os.path.exists(p) else {"pass": False, "reason": "no result (crash?)"}

    # R4: large step
    C["R4_large_step"] = {"complete": ok["big"], "pass": ok["big"],
                          "cutbacks": (R["big"]["run"].get("cutbacks") if R["big"] else None)}

    # R5: damage near 1 (stress variant)
    if PROB_STRESS is None:
        C["R5_damage_near_1"] = {"pass": None, "reason": "not run: no stress variant given"}
        return finish(res, out)
    full = complete(rs, RUNS["stress"]["cycles"], RUNS["stress"]["spc"])
    dmin = dmax = None
    if rs is not None and rs.get("fields") is not None:
        dmin, dmax = float(rs["fields"][:, 3].min()), float(rs["fields"][:, 3].max())
    in_range = (dmin is not None and TOL["d_range"][0] <= dmin and dmax <= TOL["d_range"][1]
                and dmax >= TOL["d_reached"])
    C["R5_damage_near_1"] = {"complete": full, "d_min": dmin, "d_max": dmax, "pass": bool(full and in_range),
                             "cutbacks": (rs["run"].get("cutbacks") if rs else None)}
    return finish(res, out)


def finish(res, out):
    C = res["checks"]
    res["not_run"] = sorted(k for k, v in C.items() if v.get("pass") is None)
    res["passed"] = sorted(k for k, v in C.items() if v.get("pass"))
    res["failed"] = sorted(k for k, v in C.items() if v.get("pass") is not None and not v.get("pass"))
    json.dump(res, open(os.path.join(out, "score_test1.json"), "w"), indent=1, default=float)
    write_md(res, os.path.join(out, "score_test1.md"))
    print(json.dumps({"passed": res["passed"], "failed": res["failed"], "not_run": res["not_run"]}))
    return res


def write_md(res, path):
    C = res["checks"]
    L = [f"# Test 1 score: {res['workspace']}", "",
         f"Problem {res['problem']}, {res['cycles']} cycles. Tolerances: `{json.dumps(res['tolerances'])}`", "",
         "| check | pass | detail |", "|---|---|---|"]
    for k in sorted(C):
        v = C[k]
        if "reason" in v:
            det = v["reason"]
        elif k in ("C1_key_observables", "C1_key_fields"):
            det = ", ".join(f"{q} {x:.2e}" for q, x in v["n32"].items())
        elif k == "C2_converges_to_key":
            det = ", ".join(f"{q} {x['n16']:.1e}->{x['n32']:.1e}" for q, x in v["per_quantity"].items())
        elif k == "C3_distortion":
            det = "; ".join(f"{m}: " + (", ".join(f"{q} {x[q]:.2e}" for q in OBS) if "d_max" in x else x.get("reason", ""))
                            for m, x in v["per_mesh"].items())
        elif k == "R1_history":
            det = f"declared {v['declared']}, worst change {v['worst']}"
        elif k == "R5_damage_near_1":
            det = f"complete {v['complete']}, d in [{v['d_min']}, {v['d_max']}], cutbacks {v['cutbacks']}"
        else:
            det = ", ".join(f"{a}={b}" for a, b in v.items() if a != "pass" and not isinstance(b, (dict, list)))
        L.append(f"| {k} | {'not run' if v.get('pass') is None else 'yes' if v.get('pass') else '**no**'} | {det} |")
    L += ["", "Info (not scored):", "", "```", json.dumps(res["info"], indent=1, default=float), "```"]
    open(path, "w").write("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["run", "score"])
    ap.add_argument("workspace")
    ap.add_argument("out")
    ap.add_argument("--problem", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--stress-problem", default=None)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    global PROB, PROB_STRESS, KEY
    PROB, KEY = os.path.abspath(a.problem), os.path.abspath(a.key)
    PROB_STRESS = os.path.abspath(a.stress_problem) if a.stress_problem else None
    ws, out = os.path.abspath(a.workspace), os.path.abspath(a.out)
    if a.mode == "run":
        run(ws, out, a.jobs, only=(a.only.split(",") if a.only else None))
    else:
        score(ws, out)


if __name__ == "__main__":
    main()
