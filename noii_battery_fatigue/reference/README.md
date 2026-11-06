# The reference solution (the answer key)

An independent dolfinx 0.10 implementation of the model in `../task/PUBLIC_TASK.md` §1, written by the organizer (Claude)
before the runs. None of it was shown to the agents. It follows the model's equations and 2D conventions as the task states
them, which come from the source authors' script, but contains none of their code.

- **Fields:** displacement u (quadratic), concentration c and damage d (linear), and an auxiliary field q = g(d) p_eff
  (linear), tied to g p_eff by an L2 projection so that the stress-gradient flux uses ∇q (formulation F1). F2, one
  polynomial order higher (u cubic, the others quadratic), is the check: on D0, F1 at 128 × 128 and F2 at 64 × 64 differ
  by 9.3e-5 in d_max and 8.2e-5 in ∫d.
- **Time and history:** backward Euler; ψ⁺, W, α and H at quadrature points, updated only after a converged step; the damage
  equation of a step uses H of the last accepted step.
- **Solver:** PETSc SNES Newton with LU (MUMPS).

| File | What it is |
|---|---|
| `solver.py` | The model, its weak form and the step with history update |
| `problem.py` | Problem set-up, boundary conditions, observables and the output on the host meshes |
| `d0.py` | The practice problem D0 |
| `run_d0.py`, `params_d0.json` | The D0 runs (mesh and time-step series) and their parameters |
| `mms.py` | Manufactured-solution verification (steady and transient), with results in `results/mms/` |
| `results/D0/tables.md`, `convergence.json` | D0 mesh, time-step and formulation studies |

`d0.py`, `run_d0.py` and `params_d0.json` are the organizer's files with the hidden problem's definition removed; in
`run_d0.py` the fixed location of the dolfinx environment also became the variable `FENICSX_BIN`. Nothing else in them
changed (`make_package.py` in the private repository asserts each cut). The other files are unchanged.

## Run it

```bash
FENICSX_BIN=<dolfinx env>/bin CORES=8 python run_d0.py D0:F1_u64_div40 D0:F1_u64_div80
python mms.py all
```

The two runs above give the key the scorer uses (`../scoring/reference_D0/` keeps their cycle-end values). On D0 the
64 × 64 key differs from the 128 × 128 run by 2.5e-4 in ∫d and 7.1e-4 in d_max at the last cycle end, relative to the
values (`results/D0/tables.md`).

## Checked against an agent's own solution

Run A, written without access to the key, also built its own quadratic (Q2) reference. On D0 (5 cycles at T/40) it
agrees with the key to 0.1 % in d_max and ∫d, so the damage gap of the 4-node elements (`../results/README.md`) is their discretization, not
an error in the key.
