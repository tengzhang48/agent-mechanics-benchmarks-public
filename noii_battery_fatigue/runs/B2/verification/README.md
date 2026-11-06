# Verification and reproduction

The submitted element is the four-source `enhanced_accepted_H` formulation
in [../element/uel_contract.json](../element/uel_contract.json). Its design,
time split and finite-mesh limitations are in [../design.md](../design.md).
The automatically assembled numerical tables are in [results.md](results.md),
with complete machine-readable records in [results.json](results.json).
Those tables select the exact submitted Fortran hashes. Older successful
and failed formulations are exploration evidence, not final-source results.

All commands run from the workspace root, offline, in the organizer's
prescribed Python/compiler environment. `env.sh` sets its PATH/PYTHONPATH,
single-thread BLAS and workspace-local temporary/cache directories. No
packages were installed. Host and generator package were not edited. No
Abaqus execution or independent reviewer is claimed.

## Fast compiled and analytical checks

```bash
source env.sh
python verification/check_element.py
python verification/check_accumulator.py
python verification/check_rollback.py
python verification/check_limits.py
python verification/check_history_time.py
python verification/check_energy_bridge.py
python verification/check_damaged_diffusion.py
python verification/check_mutations.py
python verification/high_order_reference.py --check
python verification/high_order_reference.py --check --qorder 6
python verification/compare_fields.py
```

- `check_element.py` compares compiled residual and state with an independent
  matrix implementation of the public laws, then finite-differences all 36
  field tangent blocks. The local modes are computed independently by
  weighted least squares. Checks include non-affine mapping, pure shear,
  biaxial tension/compression, rigid motion, compatible linear swelling,
  fatigue thresholds/floor, 24 near-failure states, noncommitting requests,
  unsupported requests, invalid data and `DTIME=0`. The most sensitive small
  tangent block has relative difference about 3.4e-6; residual agreement is
  2.3e-16. Scalar state error is below 8e-13.
- `check_accumulator.py` compiles a known-broken arithmetic ordering and
  reproduces one-ulp alpha decreases. Parenthesizing the positive energy
  increment removes them. It checks the actual monotonicity/rejection guard.
- `check_rollback.py` forces the same rejected increment before and after
  two genuine UEL trial evaluations, using the normal host stepping. Every
  saved nodal/state array agrees exactly between the two complete runs.
- `check_limits.py` checks a Fourier diffusion mode and conservation at
  beta=0 on structured/distorted meshes, time refinement, and global affine
  patches with H~1e16 and intactness below 1e-16. An archived five-field
  compiled control exhibits the original `1-d` cancellation failure.
- `check_mutations.py` deliberately changes pressure `/3` to `/2`, and
  dev_2 to an embedded dev_3. Both incorrect models pass the generator's
  derivative consistency check; the independent assembled oracle rejects
  both. These controls use the archived full-history element because their
  purpose is constitutive transcription, not the selected time split.
- `check_history_time.py` imposes a uniform linear shear ramp with beta=0
  and inactive fatigue. The exact continuum H, d and accumulated positive
  energy are analytical. Refinement from 16 to 256 steps verifies the
  one-increment damage delay and first-order convergence of d and alpha.
- `check_energy_bridge.py` uses an exactly stress-free compatible swelling
  field. The host Q1 energy equals the independently derived nonzero
  `(beta*a)^2*(3mu+lambda)*h^2/24`, verifying the reconstruction error and its
  second-order mesh dependence without relying on a simulation reference.
- `check_damaged_diffusion.py` diagnoses low-mobility concentration accuracy.
  Compiled n=8/16 fully damaged beta=0 step-charging runs match an independent
  1D reduction within 6e-15. Its n=8..2048 spatial refinement uses a banded 1D
  solve and a verified gamma-time erfc reference for continuous-space backward
  Euler. It exposes coarse consistent-mass undershoot, which resolves with
  refinement. It does not claim full 2D execution at those fine resolutions
  or prove that every severe-case discrepancy has this single cause.
- The reference check validates physical Qp Hessians, a nonzero analytical
  pressure-gradient flux, and its assembled tangent. `compare_fields.py`
  separately validates physical-space Q2/Q3/Q4 interpolation using known
  polynomials on distorted and structured meshes.

JSON files with matching check names retain actual numeric errors and
scopes. [self_review.md](self_review.md) records the single-agent review.

## Deterministic source and official entrypoint

```bash
source env.sh
python element/build.py --out tmp/rebuild
cmp element/generated_kernel.for tmp/rebuild/generated_kernel.for
cmp element/uel_contract.json tmp/rebuild/uel_contract.json
python host/run_problem.py . --mesh distorted --n 16 --seed 4 \
  --cycles 1 --steps-per-cycle 40 --out tmp/official_reproduction
```

The source is generated from handwritten Python equations after checking
four material regimes in both the full and selected declarations. The build
then inserts the explicitly documented local condensation hooks. External
handwritten sources implement condensation, validation, and the row-sum
damage reaction. [provenance.json](provenance.json) gives exact hashes,
versions, transformations and physical line counts: approximately 94% of
submitted Fortran lines originate from the generator. This count includes
unused standard helper routines; the formulation and Python equations are
handwritten. [generation.log](generation.log) and
[regeneration.json](regeneration.json) record the final byte-identical rebuild.

The official entrypoint output is preserved in
`runs/D0_enhanced_accepted_H_official_distorted_16_s4_spc40_cy1/`.
Its original host summary is `official_host_run.json`; the adjacent summary
adds provenance only. Compiler build products remain under `tmp/`.

## Global studies

```bash
source env.sh
python verification/run_case.py --n 16 --cycles 20 --spc 80 --solver-stats
python verification/run_case.py --n 32 --cycles 20 --spc 80 --solver-stats
python verification/run_case.py --n 64 --cycles 20 --spc 80
python verification/run_case.py --n 32 --cycles 20 --spc 40 --solver-stats
python verification/run_case.py --n 32 --cycles 20 --spc 160 --solver-stats
python verification/run_case.py --n 32 --cycles 20 --spc 320 --solver-stats
python verification/run_case.py --n 32 --cycles 5 --spc 320 --solver-stats
python verification/run_case.py --n 32 --mesh distorted --seed 1 --cycles 5 --spc 80
python verification/run_case.py --n 32 --cycles 5 --spc 80 \
  --problem verification/severe_problem.json --solver-stats
python verification/run_case.py --n 32 --cycles 2 --spc 160 \
  --problem verification/severe_problem.json --solver-stats
python verification/run_case.py --n 16 --mesh distorted --seed 2 --cycles 5 --spc 80 \
  --problem verification/notch_problem.json
python verification/run_case.py --n 16 --cycles 20 --spc 10 --solver-stats
```

The harness uses the unmodified host with one monolithic field group, default
tol=1e-8, line search and ten cutbacks. `--solver-stats` makes the identical
call with passive counters around the same linear solver; it does not change
tolerances, Newton, boundary conditions or stepping. Exact sources are copied
before compilation and hashed into each summary. `observables.csv` is written
at completion, fields/state at every cycle, and solver/Newton attempt records
are preserved. H and alpha monotonicity are asserted at every accepted
subincrement, not merely at saved cycles.

The severe case changes Gc to 1e-7 and alpha_T to 1e-6. The notched case
changes charging edge and displacement supports and prescribes a d=1 notch.
Neither modifies the element or embeds known solutions in it. The quarter
step-count case tests solver survival; its appreciable time error is reported.

## Continuum and independent reference comparisons

```bash
source env.sh
python verification/manufactured.py --n 8 16 32 64 \
  --mesh structured distorted --out verification/enhanced_manufactured
python verification/high_order_reference.py --n 8 --order 3 --cycles 20 --spc 80
python verification/high_order_reference.py --n 8 --order 4 --cycles 20 --spc 80
python verification/high_order_reference.py --n 8 --order 3 --cycles 20 --spc 160
python verification/high_order_reference.py --n 16 --order 3 --cycles 2 --spc 160 \
  --problem verification/severe_problem.json --lag-history
python verification/high_order_reference.py --n 16 --order 3 --cycles 2 --spc 160 \
  --problem verification/severe_problem.json --lag-history --qorder 6
```

The manufactured problem differentiates continuum stress, transport flux and
damage sources independently in SymPy. The same UEL solves them under the
normal host with an independently integrated load subtraction. The selected
element has approximately second-order displacement/damage, first-order
concentration, and decreasing complete flux error. Raw grad(p) is not strongly
convergent on the distorted family; this is a known limitation.

The NumPy high-order reference has its own constitutive/history equations,
physical shape Hessians and direct grad(g p_eff). It uses no generator code,
pressure projection, intactness alias, local enhanced modes or reaction
lumping. The generic host Newton/stepping machinery is shared. Q3/Q4 and mesh
comparisons assess reference error; it is not the unavailable organizer's
answer key. It covers top charging without a notch. `--lag-history` uses the
same accepted-history time split for severe spatial comparisons; ordinary
D0 references use current endpoint history. These differences and finite
time errors must be considered when interpreting comparisons.

Reference integration defaults to `(p+1)×(p+1)` Gauss points. Nonlinear
degradation and fatigue are not integrated exactly by a finite polynomial
rule. A separate 6×6 severe reference checks this integration sensitivity;
`--qorder 6` is also covered by the analytical kinematics/flux/area/tangent
oracle. Default reference numerics were unchanged when this option was added.

The severe consistent-reaction Q3 reference has small discrete damage
overshoots: about .0044 above one at n=8, .0018 at n=16. Its damage integral
also changes by about 1% under that refinement. The reference has measurable
spatial error; it does not establish exact pointwise bounds or exact answers.
The n=32 reference completes with overshoot .00108. Against it, UEL n=32/64
concentration L2 errors decrease from .0140 to .00770 on structured meshes
and .0293 to .0149 on distorted meshes. At n=64 the scored severe-case
energy still differs by approximately 55%/51%, respectively. That remaining
accuracy limit is not removed by the successful nonlinear solves.

Physical-space primary-field and enhanced-energy comparisons can be reproduced
with `compare_fields.py --production <run-folder> --reference <Qp-folder>
--cycle <number> --out <workspace-relative-json>`. The scored elastic energy
is always the host's Q1 reconstruction. Enhanced strain affects mechanics
and history but is invisible to that reconstruction; coarse-mesh energy
therefore has a substantial additional error that needs mesh refinement.
Comparisons additionally report reference fields sampled onto the production
nodes and evaluated with the same Q1 energy reconstruction. This diagnoses
observation discretization separately while retaining the full host error.

## Audit and remaining scope

```bash
source env.sh
python verification/audit_outputs.py
python verification/compare_enhancement.py
python verification/provenance.py
python verification/summarize.py
python verification/final_manifest.py  # after every planned source study completes
```

The output audit checks coverage, unique/connectivity indices, all saved
state finiteness, cycle history monotonicity and the intactness alias, and
independently reconstructs every saved-cycle host integral and energy. It
also checks source snapshots against recorded hashes and the host against
its earlier manifest. It excludes predecessors explicitly and reports failed
final-source runs separately. Pending jobs are shown by the summary script.

`compare_enhancement.py` checks a matched twenty-cycle n=16 study at
40 steps/cycle. The ordinary-strain and enhanced-strain source kernels are
byte-identical after removing only the enhancement hooks; guard and reaction
are byte-identical. Both use the corrected accumulator and the same properties,
history split and tolerances. `enhancement_comparison.json` retains this
control and its same-time-scheme independent Q4 comparison when available.

The accepted-history scheme has first-order temporal error. Damage bounds
follow from a discrete maximum principle on uniform square meshes; arbitrary
distortion/aspect ratios have no such guarantee. Tested distortion, loading,
and reference scopes are finite. Hidden problems, the organizer's independent
reference and actual Abaqus execution remain unverified.
