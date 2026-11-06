# Verification and reproduction

Use only the pinned environment and this workspace. No packages are installed.
The supplied `host/` is unchanged. All compilation and temporary outputs are
under `tmp/`; saved results are under this directory.

```bash
export PATH=<env>/bin:$PATH
export TMPDIR=$PWD/tmp
export PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1
mkdir -p tmp
```

The selected element has source hash beginning `ba0d3006baeb`, and uses
lagged accepted H in its damage equation. Its returned histories always use
current fields. The earlier fully implicit source, hash `fcd4c37a1a12`, is
archived in `alternatives/uel_implicit.f90`. `uel_lagged.f90` is an exact copy
of the selected source used for the controlled comparison.

## Tests

- `python verification/test_element.py`: compiled Jacobian, nonaffine analytic
  patches, rigid motion, conservation, analytical cyclic history, zero time,
  invalid orientation, repeatability after trials, and H=1e22.
- `python verification/test_reference.py`: the independent FV Jacobian,
  including the selected lagged split, and a quadratic-displacement
  counterexample to a broken Q1 pressure gradient.
- `python verification/manufactured.py --sizes 8 16 32 64`: continuum-source
  manufactured solution on structured and randomly distorted meshes.
  `--lumping 0` and `--pressure 1` reproduce spatial alternatives.
- `python verification/global_checks.py rollback`: a rejected increment before
  versus after high-history trial calls. Checks bitwise identity.
- `python verification/global_checks.py statistics`: D0, n=16, 20 cycles,
  including actual sparse linear-solve counts per Newton attempt.
- `python verification/global_checks.py large`: distorted n=16, 20 cycles,
  10 nominal steps/cycle; the ordinary host cutback mechanism is active.
- `python verification/global_checks.py saturation` and `extreme`: a constrained
  notch, then the same problem with a very small fatigue threshold, to test
  near-one and fully saturated damage. These are survival tests; temporal
  accuracy in the extreme threshold case is not established.
- `python verification/uniform_history.py`: homogeneous fatigue/damage, with an
  independent scalar-root implicit reference at 5120 steps/cycle. It exposes
  coarse-time history error without spatial error.
- `python verification/q2_reference.py --test` and `--n 8` (also n=4,16,32):
  a separate higher-order 2D formulation with direct Q2 displacement Hessians,
  no pressure/source auxiliary fields, and 5x5 quadrature. It uses the generic
  unchanged host with reference-only nine-node macroelements; the scored
  source/contract remain four-node. `--n 16 --spc 80` and `160` refine time.

The custom checks use the same unmodified one-block Newton/line search and
transactional host. They neither stagger fields nor increase its normal
25-iteration/10-cutback limits. The fixed body sources of the manufactured
test are subtracted in a test-only batch adapter. The selected UEL has no
test-dependent code.

## Full-domain cyclic studies

For example:

```bash
python verification/run_host.py selected_D0_s16_spc40 --mesh structured --n 16 --cycles 20 --steps-per-cycle 40
python verification/run_host.py selected_D0_d16_spc40 --mesh distorted --seed 1 --n 16 --cycles 20 --steps-per-cycle 40
```

`run_host.py` invokes the original `host/run_problem.py`, builds under
`tmp/runs/NAME/build`, and copies its numerical results to `verification/runs`.
Each new run includes an input contract/source-hash manifest. Mesh sizes are
8,16,32 (20 cycles), with structured 64 (five cycles). Time studies use
40,80,160 steps/cycle on n=16 structured and distorted meshes. All other
parameters are the original D0 file.

The directories prefixed `selected_` contain the delivered formulation.
Earlier directories prefixed `final_` retain their historical names but
belong to the old **implicit** candidate. Failed attempts are preserved in
`global_checks/`, with source hash, full problem, attempts and accepted commits.

## Independent coupled reference

`reference_1d.py` eliminates displacement analytically in a one-dimensional
constrained strip and uses cell-centred finite volumes for c,d. It does not
call the UEL for reference physics and uses no q,z variables.

```bash
python verification/reference_1d.py fv --lag --n 512 --spc 80
python verification/reference_1d.py uel --n 128 --spc 80
python verification/reference_1d.py fv --n 256 --spc 320
```

The `fv` form is fully implicit in current history. `fv --lag` uses the same
first-order history split as the selected UEL, so spatial error can be
isolated. Both reference time integrations approach the same PDE under dt
refinement. `uel` solves a 2D strip with the host; its integrals are divided by
strip width. Directories `implicit_uel_*` retain the old UEL comparison.

## Read results

- [RESULTS.md](RESULTS.md): generated tables and limits.
- `results.json`: the same evidence in structured form.
- `runs/*/observables.csv`: all accepted nominal output times, including the
  host's threshold-sensitive damage area.
- `field_convergence.json`: common-physical-point Q1 field comparisons against
  n=64 (a numerical refinement study, not an independent 2D solution).
- `q2_comparison.json`: delivered Q1 versus the separately discretized Q2
  family, including Q2 self-convergence and actual flux comparisons.
- `provenance.json`: environment and source/host hashes.

Once all cyclic runs are present, regenerate field comparisons and tables:

```bash
python verification/field_convergence.py
python verification/compare_q2.py
python verification/summarize.py
```

The pressure gradient alone remains noisy on random distortions; full flux
and primary fields converge in the manufactured test. The first-order
history lag requires time refinement, especially near abrupt fatigue growth.
Hidden reference problems and actual Abaqus execution have not been tested.
