# Research log

All times are UTC from the system `date` command. Single agent; no internet,
installs, outside-workspace challenge reads, or Abaqus execution.

- 2026-10-20T23:10:09Z — Initial clock reading after reading ORGANIZER_NOTE.md,
  PUBLIC_TASK.md and host/README.md. Arm B, pinned abaqus_ufl a02ed7c.
  Read the bundled AI skill, API guide and example verification instructions.
  `git status --short` reports this workspace is not a Git worktree.
- Initial formulation decision: equal-order isoparametric Q1 fields u, c, d and
  an additional scalar p, constrained by a global L2 projection of g(d)p_eff.
  This makes the transport flux meaningful without elementwise displacement
  Hessians. Plan mesh/distortion and time refinement to assess the projection.
- Planned history scheme: implicit trial alpha and H recomputed from accepted
  H, alpha, W at every quadrature point. Fully differentiate active branches.
  Return new state to the host, which alone accepts or rolls it back.
- Tool scope: use the standard Quad4 UEL generator, Python constitutive
  declarations and generated complex-step tangents. A small handwritten
  standard-signature wrapper will add validation and a DTIME=0 rejection.
  No changes planned to host/ or the abaqus_ufl package.
- Review will be performed by this single agent, as organizer conditions
  override the skill's suggestion of independent agents.
- 2026-10-20T23:13:23Z — Generated 5-DOF prototype passed compiled, independently
  assembled residual/history and finite-difference tangent checks. Pure-shear
  analytic oracle agrees to 3.6e-15. Initial generator failure: `dict(...)` in a
  state return is not translated; replaced by a literal dictionary, with no
  package change. D0 structured n=8 and n=16, 2 cycles x 40 steps, completed.
  n=8 showed a nodal damage overshoot to 1.063; n=16 reached 0.291 with no
  cutbacks. Retain this as evidence of underresolution, not convergence.
- 2026-10-20T23:16:01Z — Exhausted-fatigue scalar Newton test (alpha_old=.02,
  pure shear .01) exposed cancellation in 1-d at H~1.3e16: residual stalled
  near 1e-1 after d rounded to 1. Pivot: introduce intactness s=1-d as a solved
  field with a mass constraint, put its stable Helmholtz equation on the d
  rows, and use s in degradation. This is algebraically identical in exact
  arithmetic; new final contract is 6 DOFs per node. The 5-field code is
  retained only in tmp/initial5; its n=16 result is labeled
  verification/exploration/initial5_*.
- 2026-10-20T23:17:49Z — Six-field intactness element passes the independent
  assembled checks. Exhausted-fatigue Newton converges in two updates, with
  residual 1.44e-16 and intactness 7.47e-17. Manufactured forcing test initially
  failed: initializing the extra intactness to zero while applying a full
  independent body force gave a kappa-stiff initial tangent and an invalid
  displacement trial. Set the verification initial state to its physical
  intactness s=1. D0 retains the host's zero extra-DOF initialization and works.
  The original failure trace is in verification/exploration/.
- 2026-10-20T23:21:45Z — Independent smooth manufactured solution, structured
  and distorted n=8,16,32,64. Structured c L2 error falls from 1.00e-3 to
  1.51e-5; distorted c from 1.75e-3 to 2.23e-4. Distorted raw grad p error
  stays near 0.033, while the actual flux error falls from 2.17e-3 to 5.25e-4.
  Do not claim strong pressure-gradient convergence. Exact rollback equality
  succeeds on a distorted mesh, including all nominal output states.
- 2026-10-20T23:22:35Z — Tightened mechanics and pressure residual scales from
  1 to .01 to retain equilibrium accuracy at g near kappa. Added explicit
  rejection for unsupported LFLAGS and negative/nonfinite history. Repeated
  compiled and forced-rollback checks. Final kernel equations unchanged.
- 2026-10-20T23:24:47Z — Five-cycle low-toughness run at n=16 completes, with
  two rejected attempts and nodal d overshoot up to 1.015. A distorted notched
  run at n=16 completes five cycles without cutbacks. Global Fourier diffusion
  matches the exact Q1 discrete eigenmode to 2.8e-16; mass conserved to
  1.1e-16. Fully damaged global affine patches on distorted n=8,16,32 give
  displacement error below 1.5e-17 at H=1.11e16.
- 2026-10-20T23:25:52Z — Explored a fully mixed hydrostatic stress alternative,
  declared in verification/exploration/mixed_stress.py. It puts 1.5 p I into
  momentum and uses the compliance-weighted constraint p/(g C)-tr(eps_e)=0.
  Manufactured concentration error improves modestly but remains first order
  on distorted meshes and grad p still lacks strong convergence. Retain the
  simpler primary element; do not imply this alternative was globally tested.
- 2026-10-20T23:27:51Z — Built an independent NumPy high-order Qp reference,
  with physical shape Hessians and the strong stress gradient, no abaqus_ufl
  import. Quadratic physical kinematics on a non-affine cell agree to 5.0e-16;
  analytical nonzero pressure-gradient flux agrees to 5.5e-17; assembled
  tangent relative error 9.5e-11. Q2 n=8 D0 completes 2 cycles x 40 steps,
  final d_max=.190601 and D_int=.0668976. More reference refinement is running.
- 2026-10-20T23:32Z–23:47Z — Global low-toughness n=32 studies exposed
  repeated Newton failures in the fully implicit history formulation. Saved
  failed trial states, solver reports and kernel snapshots under
  `verification/exploration/`. Tested nodal/corner history lumping, nodal
  reaction lumping with Gauss history, a smaller complex-step perturbation,
  a blended history-switch tangent, phase residual rescaling, sqrt(g) row
  scaling, and removal of the irrelevant det(I+grad u) guard. None resolved
  the severe n=32 failure with current-iterate history; first failures occurred
  around steps 12–17 after ten cutbacks. Central differences crossed active
  history branches, and line searches found no residual-decreasing length.
  Extended-precision interpolation checks ruled out cancellation in the
  Gauss intactness evaluation as the cause at those trial states.
- 2026-10-20T23:47Z–23:59Z — Pivoted to start-of-increment H in the damage
  equation, with exact current-field endpoint H/alpha/W returned to the host.
  This is an explicitly documented first-order history time split; fields
  still solve together in the host's single monolithic Newton group. Added
  row-sum reaction lumping to improve damage bounds. This combination completed
  the severe n=32 case without cutbacks. Large-step D0 also completed, but
  accuracy is visibly time-step dependent. Launched long mesh/time studies.
  These basic Q1 results showed excessive cyclic damage on coarse meshes:
  n=16 reached d_max=.953 after twenty cycles, versus approximately .294 for
  the independent strong-gradient Q3 reference. Solver success alone is
  insufficient evidence of accuracy.
- 2026-10-20T23:59Z — Independent Q3/Q4 reference refinement and twenty-cycle
  runs in progress. The reference uses direct physical shape Hessians and
  the public constitutive equations, no generator, projected pressure or
  intactness field. It is a validation implementation, not an answer key.
- 2026-10-21T00:06Z–00:14Z — Investigated incompatible/enhanced strain modes
  to reduce parasitic Q1 swelling stress. Added four element-local quadratic
  displacement modes with zero volume-average physical gradients and solved
  their elastic equilibrium exactly. Derived and implemented the full
  condensed Schur correction for mechanics, concentration and pressure
  couplings. No local mode or factorization persists between calls. An
  independently assembled residual/tangent and non-affine affine-swelling
  patches passed. A hybrid feedback of projected pressure into stress did not
  improve manufactured concentration convergence and was not selected.
  Enhanced strain dramatically improved D0 and severe-load agreement with
  the independent reference; promoted it after further compiled checks.
- 2026-10-21T00:14Z–00:16Z — A 320-step/cycle severe run unexpectedly failed
  its history monotonicity guard. Isolated a genuine floating-point bug:
  `alpha_old + W_new - W_old` can be one ulp below alpha_old despite
  W_new > W_old. Changed this to `alpha_old + (W_new - W_old)`. Constructed
  numerical witnesses that reject the predecessor and pass the correction;
  the exact same formerly failing basic run now completes 640 increments
  with no cutback. The fix is arithmetic grouping, not a history clamp.
- 2026-10-21T00:23Z — Froze the enhanced, accepted-history source with the
  accumulator correction. Repeated element, rollback, Fourier, conservation,
  and fully damaged affine checks on this source. Fresh authoritative global
  runs use the `enhanced_accepted_H` formulation tag and preserve actual
  source bytes before compilation. Earlier variants remain exploration
  evidence and are not claimed as final-source runs.
- 2026-10-21T00:28Z — Final source completed twenty D0 cycles, structured
  n=16, 80 steps/cycle: d_max=.285369, D_int=.0724268, C_int=.512520,
  1600 accepted increments, no rejection, exactly two Newton updates each.
  Rewrote design.md around the selected formulation and its limitations.
- 2026-10-21T00:33:53Z — Clock check. All final manufactured n=8,16,32,64
  structured/distorted runs completed. Displacement/damage are approximately
  second order; concentration is approximately first order. On distorted
  meshes raw grad(p) error stays near .034, while complete flux error falls
  from .00217 to .000438. Record this weakness explicitly. Final large-step
  twenty-cycle and distorted notched five-cycle tests also completed without
  cutback. Longer n=32/n=64 and independent-reference studies continue.
- 2026-10-21T00:35:19Z — Final kernel and contract regenerate byte-identically.
  The unmodified official host entrypoint completes distorted n=16, seed=4,
  one cycle x 40 steps without cutback. Preserve its original summary and
  fields, alongside source provenance. Known-wrong generated pressure and
  deviator controls are again rejected by the independent assembled oracle.
- 2026-10-21T00:39:24Z — Final severe n=32 completes five cycles x 80 steps,
  d_max=.999882, D_int=.368691, no rejected increment, 2–3 Newton updates.
  Q3/Q4 independent twenty-cycle references differ by 2.35% in final maximum
  damage and .55% in damage integral, so reference error is not negligible.
  Launch further Q3 mesh and Q4 time refinement instead of treating either
  as an exact solution. Added independent physical-space field comparison;
  its polynomial oracle first exposed an inverse-map assertion on wrong
  distorted-cell candidates. Corrected candidate rejection, then verified
  physical quadratic reproduction to below 2e-15 for Q2/Q3/Q4.
- 2026-10-21T00:41Z — Analytical uniform-shear oracle completed on distorted
  n=8, 16/32/64/128/256 time steps, with inactive fatigue. It gives exact
  H(t)=10t^2, an independent continuum damage solution and an analytical
  maximum of W for accumulated alpha. Damage and alpha errors approach order
  one (finest rates 1.0075 and 1.0133). The discrete one-increment damage
  delay matches the analytical split formula to 3.4e-16. All accepted
  histories remain monotone. Added reproducibility and single-agent-review
  documents; long mesh/time/reference runs continue.
- 2026-10-21T00:45:58Z — Derived the host Q1 energy for an exactly compatible
  stress-free swelling field: `(beta*a)^2*(3mu+lambda)*h^2/24` on the unit
  square. Direct host evaluation at n=8,16,32,64 matches to relative error
  below 8e-15. This verifies a real O(h^2) reconstruction error independently
  of the reference solver. In a severe n=32 comparison, the reference sampled
  on production nodes yields Q1 energy 2.239e-8 versus the UEL's 2.194e-8,
  while the resolved reference energy is 5.503e-9. Do not subtract this
  observation error or hide the large coarse-mesh scored-energy discrepancy.
- 2026-10-21T00:47:59Z — Final severe n=32, two cycles x 160 steps, completes
  without cutback. At the same dt, independent Q3 n=16 differences are .56%
  in D_int and -.030% in C_int; primary damage L2 error is .0137. A previous
  comparison at differing dt had much larger field discrepancies, motivating
  an additional independent severe time refinement. Started distorted severe
  n=32 as a separate robustness check.
- 2026-10-21T00:49:52Z — Final n=32 D0, twenty cycles x 40 steps, completes
  without cutback, d_max=.268185. Final n=32 at 80 steps/cycle also completes
  shortly afterward, d_max=.287363, D_int=.0711111, C_int=.512480. Histories
  remain monotone at every accepted subincrement. Auditing those final-source
  outputs reconstructs all saved-cycle observables to roundoff. Launch a
  twenty-cycle 320-step run to extend the global time-refinement sequence.
- 2026-10-21T01:01:48Z — A controlled late alternative combines current-H Gauss damage
  with the same enhanced strain. Its independent non-affine residual/state
  and full condensed tangent pass (normalized FD error 2.34e-9). An initial
  invalid FD harness reused the private driver's in-place SVARS buffer;
  copying incoming state for every perturbation corrected that test without
  a Fortran change. D0 n=16, two cycles x 40 steps, completes (d_max=.19609),
  but severe n=32 fails at step 19 after ten cutbacks, d_max=1.02530, phase
  residual 5.60e-7. This directly tests whether enhancement alone fixes the
  earlier implicit-history issue; it does not. Keep the frozen accepted-H
  lumped formulation. Experimental source and failure evidence are preserved.
- 2026-10-21T01:12:00Z — Matched enhancement control completed: structured n=16,
  twenty cycles x 40 steps, identical corrected accumulator, history split,
  reaction, parameters and residual scales. Removing only the enhancement
  hooks makes generated kernels byte-identical; guard/reaction are already
  byte-identical. Ordinary Q1 strain reaches d_max=.953243, D_int=.135451;
  enhanced strain reaches .267728, D_int=.0700064. Both complete without
  cutback. Added an independent Q4 reference at the same accepted-H time
  scheme and dt to assess this controlled spatial comparison.
- 2026-10-21T01:28:27Z — The same-split independent Q4 twenty-cycle, 40-step
  reference completed without cutback: d_max=.267448, D_int=.0682562.
  Matched n=16 enhanced d_max differs by .105%, D_int by 2.56%; ordinary
  Q1 d_max is 256% high. Host Q1 enhanced energy remains 50.1% high; peak
  agreement must not be presented as full-field accuracy. The expanded
  output audit independently reconstructs all five saved-cycle observables
  (d_max, area above .5, D_int, C_int, E_el) within CSV precision.
  Twenty-cycle n=32 time studies at 40/80/160 steps are complete without
  cutback (d_max=.268185/.287363/.298418); 320 steps continues. Severe n=64
  completed two cycles x 160 steps through d_max=.999941 without cutback.
  A concentration-error plateau relative to the older severe reference
  motivates one further Q3 mesh and 6x6 quadrature refinement of that
  reference. The new optional quadrature setting passes its analytical
  geometry/flux/area/tangent oracle; default reference numerics are unchanged.
- 2026-10-21T02:00:27Z — Diagnosed a separate transport resolution limitation.
  With beta=0 and uniform accepted H=1e16, the element reduces to diffusion
  with kappa*D=1.4e-6. Compiled n=8/16 step-charging solutions match an
  independent 1D consistent-mass BE solve within 5.8e-15. At t=8 the
  diffusion length is .003347, far below h at n=64. The coarse 1D profile
  undershoots c0=.35 (n=64 final minimum .280658, L2 error .01394); fine
  1D n=2048 resolves it (L2 1.56e-5, no undershoot). This fine refinement
  is not a full 2D host run. A gamma-time erfc formula provides the exact
  continuous-space BE reference: its resolvent recurrence and 64-vs-96
  integration check agree within 7e-16. Report concentration positivity
  and low-mobility field resolution as limitations; no clipping or artificial
  diffusivity was added. This mechanism is not a complete attribution of
  all coupled severe-case discrepancies.
- 2026-10-21T02:39:55Z — Twenty-cycle n=32 time sequence is complete: 40/80/160/320
  steps per cycle give d_max=.268185/.287363/.298418/.304703. All four
  runs have zero cutbacks and nondecreasing H/alpha at every accepted step.
  Distorted n=16/32/64 five-cycle comparisons against the same accepted-H
  Q4 reference show c L2 errors .000684/.000213/.0000778 and host-energy
  relative differences 63.7%/16.3%/4.2%. This separates spatial convergence
  from the finite-step bias of endpoint-history comparisons. Q3 n=16 and
  Q4 n=8 endpoint-history twenty-cycle references agree within .76% in
  d_max, .17% in D_int and .13% in energy; neither is treated as exact.
  The sole remaining source run is n=64, twenty cycles x 80 steps. The
  finer severe reference has completed its first cycle and continues.
- 2026-10-21T02:54:37Z — The finer severe Q3 n=32 reference completes two cycles x
  160 steps without cutback. Its d_max=1.001081, D_int=.368328,
  C_int=.492795 and energy=5.37035e-9 remain finite-resolution evidence.
  Against this reference, production n=32/64 c L2 errors are .0140/.00770
  (structured) and .0293/.0149 (distorted); d L2 also decreases. The
  concentration plateau seen against the older reference is reduced, but
  neither the reference nor these mesh resolutions are exact. At n=64,
  scored severe energy is still about 55%/51% high; retain that explicit
  accuracy limitation. Refreshed every legacy field comparison to include
  the same sampled-Q1 energy diagnostic, and made unequal time steps visible
  in the report. All 18 independent reference runs are now complete.
- 2026-10-21T03:31:34Z — The final n=64 D0 study completes twenty cycles x 80 steps,
  d_max=.2856718183, with 1600 accepted increments and zero cutbacks.
  The full final-source cohort is now 21 runs and 21,960 accepted nominal
  outputs, all without cutback; H/alpha remain nondecreasing at every
  accepted increment. A fresh evidence update, after joining the old worker,
  regenerates the remaining comparisons with current diagnostic code.
  All 31 planned comparison pairs are present (35 unique comparisons total).
  Against the same-split Q4 reference, structured n=16/32/64 concentration
  L2 errors are .000230/.0000651/.0000223 and scored energy differences
  50.6%/12.8%/3.16%. Final n=64 damage L2 error is .000205. These improve
  under mesh refinement; the separate severe-case energy/transport limits,
  finite reference error, and first-order time error remain documented.
- 2026-10-21T03:32:08Z — Final completeness/integrity gate passes: 21 source-matching
  global runs, 21,960 accepted outputs, zero rejected production attempts,
  18 completed independent references, 35 unique comparisons covering all
  31 planned pairs, and 11 passed check records. The output audit checks
  4,496,916 saved scalar values. Source, contract and verification Python
  hashes match final provenance; host manifest matches its earlier record.
  Required deliverables are complete. Stop within the six-hour budget,
  retaining the documented pressure/time/concentration/energy limitations;
  hidden problems, organizer answers and Abaqus remain unverified.
