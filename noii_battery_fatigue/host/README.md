# The host: runs a submitted element on a problem

The host plays the role of Abaqus/Standard. It compiles the element named in a run's `element/uel_contract.json`, calls it
through the standard UEL argument list, and solves the global problem: one monolithic Newton iteration per increment with a
line search, start-of-increment SVARS at every call (kept only after convergence, restored on rejection), and cutbacks when
the element returns PNEWDT < 1 or Newton fails. The call convention is stated in full in `../task/host_README_as_given.md`
§1; the host keeps it exactly, including TIME(1) = TIME(2) = the time at the end of the increment and KSTEP = KINC = 1.

| File | What it does |
|---|---|
| `fehost.py` | Reads the contract and builds and calls the element through the public abaqus_ufl.fe (`build_external_uel`, `CompiledAbaqusElement.from_module`, `CompiledElementState`, `solve_steps_staggered`). |
| `run_problem.py`, `meshes.py`, `observables.py` | Run a problem file on structured or distorted 4-node meshes and write `observables.csv`, `fields_cycle_XXX.npz` and `run.json`. Unchanged from the host the agents had. |
| `dev_problem/d0.json` | The practice problem D0: the unit square without a crack, charged cyclically on its top edge for 20 cycles. |

## Setup

Python 3.10 or newer (checked with 3.13) and gfortran (the element is compiled with `numpy.f2py`), and the public
abaqus_ufl with its finite-element extras (NumPy, SciPy, Meson, Ninja) at the commit these files were checked with:

```bash
pip install "abaqus_ufl[fe] @ git+https://github.com/tengzhang48/abaqus_ufl@691c4eb"
```

## Run an element

From this folder, with a run folder that holds `element/` (for example `../runs/B2`):

```bash
python run_problem.py ../runs/B2 --mesh structured --n 16 --out out/B2_s16
python run_problem.py ../runs/B2 --mesh distorted --n 32 --seed 1 --steps-per-cycle 40 --cycles 2 --out out/B2_d32
```

## The same calculation as the host the agents had

The agents' host vendored the solver part of abaqus_ufl.fe at a development revision (commit 9ac03e8 of 27 September
2026; frozen for the challenge on 20 October 2026). Its solver
functions were ported to the public abaqus_ufl (pull requests 9 and 10, merged on 6 November 2026; checked at `691c4eb`); the call
convention of this task stays in `fehost.py`. The two hosts were compared on run B2's element on the practice problem, in
three runs that cover plain stepping, Newton failures with cutbacks, and a distorted mesh (structured n = 16 at T/40 for 4
cycles; structured n = 16 at T/10 for 4 cycles with the inner Newton limited to 2 iterations, which forces 4 cutbacks;
distorted n = 16, seed 1, at T/40 for 2 cycles). All 88 checks were exact: the same attempts and cutbacks, the same TIME,
DTIME and KINC at every element call, the same PNEWDT, bit-identical U and SVARS after every accepted increment, and equal
observables. The official scores were also recomputed on this host for run B2 and came out identical (see
`../results/README.md`).
