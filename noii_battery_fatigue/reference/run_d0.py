"""Run script for the practice problem D0 (parameters in params_d0.json).

  CORES=18 FENICSX_BIN=<fenicsx-env>/bin python run_d0.py all            # every run below, scheduled on <= CORES cores
  python run_d0.py D0:F1_u64_div40 ...     # selected runs

Each run writes results/D0/NAME/ (observables.csv, fields_grid101.npz, host_nodes.npz,
meta.json, newton_log.json, history_end.npz, fields_dofs.npz for the exported runs).
"""
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ENV = os.environ.get("FENICSX_BIN", os.path.dirname(sys.executable))  # a dolfinx 0.10 environment
CORES = int(os.environ.get("CORES", "18"))
PP = json.load(open(os.path.join(HERE, "params_d0.json")))

# (problem, name): (form, N, div, nproc, extra params, export dofs); listed in start order
RUNS = {}
for prob in ("D0",):
    RUNS.update({
        (prob, "F1_u128_div40"): ("F1", 128, 40, 4, {}, True),
        (prob, "F2_u64_div40"): ("F2", 64, 40, 2, {}, True),
        (prob, "F1_u64_div160"): ("F1", 64, 160, 2, {}, False),
        (prob, "F1_u64_div80"): ("F1", 64, 80, 1, {}, False),
        (prob, "F1_u64_div40"): ("F1", 64, 40, 1, {}, True),
        (prob, "F1_u64_div40_chi0"): ("F1", 64, 40, 1, {"chi": 0.0}, True),
        (prob, "F1_u32_div40"): ("F1", 32, 40, 1, {}, False),
        (prob, "F1_u16_div40"): ("F1", 16, 40, 1, {}, False),
        (prob, "F2_u32_div40"): ("F2", 32, 40, 1, {}, False),
        (prob, "F2_u16_div40"): ("F2", 16, 40, 1, {}, False),
        (prob, "F1_u128_div40_chi0"): ("F1", 128, 40, 4, {"chi": 0.0}, False),
    })


def cmd_for(key):
    prob, name = key
    form, n, div, npr, extra, export = RUNS[key]
    prm = dict(PP[prob])
    prm.update(extra)
    out = os.path.join(HERE, "results", prob, name)
    os.makedirs(out, exist_ok=True)
    cmd = [f"{ENV}/python", os.path.join(HERE, "d0.py"), "--problem", prob, "--n", str(n), "--form", form,
           "--div", str(div), "--ncyc", str(prm.pop("ncyc")), "--out", out, "--save_hist"]
    cmd += sum([[f"--{k}", str(v)] for k, v in prm.items()], [])
    if export:
        cmd.append("--export_dofs")
    if npr > 1:
        cmd = [f"{ENV}/mpiexec", "-n", str(npr)] + cmd
    return cmd, npr, out


def main():
    if sys.argv[1] == "all":
        keys = list(RUNS)
    else:
        keys = [tuple(a.split(":")) for a in sys.argv[1:]]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", TMPDIR=os.path.join(HERE, "tmp"))
    queue = list(keys)
    running = []
    while queue or running:
        used = sum(r[1] for r in running)
        for key in list(queue):
            npr = RUNS[key][3]
            if used + npr <= CORES or not running:
                cmd, npr, out = cmd_for(key)
                log = open(os.path.join(out, "run.log"), "w")
                pr = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=HERE)
                running.append((pr, npr, key, time.time()))
                queue.remove(key)
                used += npr
                print(time.strftime("%H:%M:%S"), "started", key, npr, "ranks", flush=True)
        time.sleep(5)
        for r in list(running):
            if r[0].poll() is not None:
                running.remove(r)
                print(time.strftime("%H:%M:%S"), "finished", r[2], "rc", r[0].returncode,
                      f"{time.time() - r[3]:.0f}s", flush=True)


if __name__ == "__main__":
    main()
