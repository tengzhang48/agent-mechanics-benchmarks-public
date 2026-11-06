# Status

Implemented and compiled the six-field Quad4 UEL in `element/`, with four
statically condensed strain modes, an auxiliary pressure projection, stable
intactness, a lumped damage reaction and first-order accepted-history time
splitting. The bundled generator produces the constitutive/history code,
bulk assembly and tangents; local condensation and the lumped reaction are
handwritten. Approximately 94% of submitted Fortran physical lines originate
from the generator, including unused standard helpers; the Python equations
and formulation are handwritten. Host and generator package are unchanged.
Work is offline and inside this workspace, using the prescribed environment.
No Abaqus run was made.

All 21 final-source global studies are complete: 21,960 accepted output
steps and zero cutbacks. All 18 independent Qp reference studies are also
complete. Exact run records, source hashes and numerical tables are in
[verification/results.md](verification/results.md).

Validated:

- Eleven compiled, analytical, transactional and integrity check records pass.
  These include independent residual/history assembly, all tangent field
  blocks, literal pressure/energy conventions, near-failure states,
  DTIME=0/invalid-state guards and noncommitting requests.
- Rejection before and after genuine UEL trial calls gives exactly identical
  complete output/state arrays. H and alpha never decrease in the completed
  final-source runs, checked at every accepted increment.
- Twenty D0 cycles complete on n=16/32/64 at 80 steps/cycle. The n=32 twenty-cycle
  time sequence at 40/80/160/320 steps completes, with final d_max
  .268185/.287363/.298418/.304703. Instrumented studies use 1–3 Newton
  updates per increment with the host's default tolerance and line search.
- Severe loading completes five cycles on n=32, and two cycles on n=32/64
  structured and distorted meshes, through d_max=.999966. A distorted
  notched problem with different supports/charging and the quarter-step-count
  twenty-cycle D0 run also complete.
- Smooth full-PDE manufactured refinement through n=64 shows approximately
  second-order displacement/damage and first-order concentration, including
  distortion. Complete transport flux error decreases.
- Against the independent Q4 reference using the same history split and time
  step, distorted n=16/32/64 concentration L2 errors are
  .000684/.000213/.0000778. Host-energy differences decrease from
  63.7% to 16.3% to 4.2%. The reference retains finite discretization error.
- The corresponding structured twenty-cycle concentration L2 errors are
  .000230/.0000651/.0000223 for n=16/32/64; host-energy differences are
  50.6%/12.8%/3.16%. At n=64, damage L2 error is .000205 and final d_max
  differs by -.285% from that finite reference.
- Against the finer severe Q3 n=32 reference, production n=32/64 concentration
  L2 errors decrease from .0140 to .00770 (structured) and .0293 to .0149
  (distorted). Damage L2 errors also decrease. This is measured refinement,
  not an exact error bound.
- Analytical Fourier diffusion and mass conservation, distorted fully damaged
  affine patches at H~1e16, first-order history-time convergence from 16 to
  256 steps, and the host's Q1 energy reconstruction formula pass.
- The unmodified official host entrypoint completes a distorted n=16 run.
  Regeneration gives byte-identical kernel/contract. Independent output audits
  reconstruct all saved-cycle observables and check source snapshots.
- A tiny-increment regression reproduces a one-ulp alpha decrease in the
  predecessor and verifies the corrected parenthesized positive increment.
  The formerly failing 320-step severe run completes after that correction.

Known limitations:

- Raw L2-projected pressure gradients do not converge strongly on the tested
  randomly distorted manufactured family. Primary concentration converges
  about first order; complete flux error decreases in the measured studies.
- Accepted-history splitting has first-order time error. Large nominal steps
  can complete with inaccurate fatigue/damage histories; time refinement is
  necessary.
- Near complete damage, diffusivity reaches kappa*D and thin concentration
  layers can remain unresolved even at n=64. A compiled beta=0 limit agrees
  with an independent 1D reduction, which exhibits consistent-mass undershoot
  on coarse meshes and converges under refinement. Concentration positivity
  is not guaranteed. Compiled checks use n=8/16; finer evidence there is 1D.
- The host reconstructs energy from Q1 displacement, while mechanics/history
  use enhanced strain. Severe-case scored energy remains about 55% high on
  structured n=64 and 51% high on distorted n=64 relative to the finite Q3
  reference. Successful nonlinear solves do not resolve this accuracy limit.
- Only uniform-square-mesh damage bounds follow from a discrete maximum
  principle. Distorted-cell bounds have a measured validation envelope.
- The independent reference has finite spatial, integration and time error;
  severe consistent-reaction Q3 can slightly exceed d=1. A separately
  generated current-H enhanced Gauss alternative fails the severe n=32
  solve at step 19 despite a verified residual/tangent. It remains a rejected
  alternative, separate from the submitted source.

Hidden problems, the organizer's reference, resolutions beyond the measured
mesh/time envelope, and production Abaqus execution remain unverified.
[design.md](design.md) and [verification/README.md](verification/README.md)
describe the formulation, rejected alternatives, reproducible commands and
scope of each check. [verification/final_manifest.json](verification/final_manifest.json)
records the final completeness and source-integrity gate.
