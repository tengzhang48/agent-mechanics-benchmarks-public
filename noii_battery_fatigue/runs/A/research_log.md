# Research log

All work is offline, confined to this workspace, in arm A. The host is used
unchanged. Sources are written directly from PUBLIC_TASK.md; no generators or
external examples are used. Compiler/Python environment:
`<env>`. Builds and other temporary files
go under `tmp/`; Python bytecode writing is disabled.

- 2026-10-20T19:38:48Z (`date -u`): Read ORGANIZER_NOTE.md first, then the
  public task and host README. Inspected the host's assembly, transactional
  state, compilation and solver conventions. The supplied D0 has 20 cycles
  and a very small fracture toughness; testing only early elastic increments
  would be insufficient.
- Initial formulation: Q1 displacement/concentration/damage with continuous
  Q1 auxiliary pressure `q = g p_eff` in an L2 constraint. This avoids broken
  second derivatives of displacement in the concentration flux. Plan to
  compare against an independent 1D reduction and manufactured patch tests,
  and study structured/distorted mesh convergence.
- Anticipated saturation issue: the fatigue floor permits H far beyond
  1/epsilon. Use an auxiliary damage source `z = (1-kappa)H(1-d)` with its
  algebraic equation divided by `1+(1-kappa)H`. This leaves the continuous
  model unchanged and avoids subtracting huge reaction terms near d=1.
  Its finite-element approximation must also be checked by refinement.
- 2026-10-20T19:41:33Z (`date -u`): First compiled-element tests passed:
  active-history tangent relative error about 2e-11, analytic constant-strain
  forces on a nonaffine quad, conservation, repeated trial calls, zero time,
  invalid orientation, and artificial H=1e22 saturation.
- 2026-10-20T19:43:27Z (`date -u`): Initial fully integrated source split
  completed D0, but coarse n=8 damage exceeded one (20-cycle maximum 1.154).
  Pivoted to nodal quadrature of d,z in their zero-order terms. H and its
  sensitivities remain at 3x3 Gauss points. The coarse maximum became 0.5283.
  The original option remains reproducible with JPROPS(1)=0.
- 2026-10-20T19:45:30Z (`date -u`): Manufactured forcing derived from the
  continuum equations using Sympy. First invocation failed because the
  staggered host API has no external_fn parameter; corrected the test to
  subtract fixed element body loads in its batch adapter, without changing
  the host. Structured fields converge approximately at order two. Randomly
  distorted meshes expose pressure-gradient noise and slower concentration
  convergence. Tested mixed volumetric stress (JPROPS(2)=1): it does not
  remove that noise; retained direct stress (default 0).
- 2026-10-20T19:47:23Z (`date -u`): Started final-formulation D0 spatial and
  temporal studies. Builds are under tmp; immutable host CLI results are
  copied into verification/runs. New runs include a source hash and contract.
- 2026-10-20T19:49:50Z (`date -u`): Global forced rejection before/after
  deliberately high-history UEL trial calls gives bitwise identical U and
  SVARS. An extreme alpha_T=1e-15 test with notch and four constrained corners
  exhausted cutbacks in the first nominal step. This is a failed attempt,
  not evidence of successful saturated-fatigue dynamics.
- 2026-10-20T19:52:03Z (`date -u`): Quarter-resolution large-step D0
  (10 steps/cycle, distorted n=16, 20 cycles) completes with 25 cutbacks.
  A constrained notch test with G_c=1e-9 also stalls. Independent 1D
  cell-centred FV reference and UEL strip solutions agree under refinement;
  the reference eliminates displacement analytically, uses no q or z, and
  reaches alpha above alpha_T.
- 2026-10-20T19:53:03Z (`date -u`): Toughness sweep without a notch:
  distorted n=8, five cycles completes for G_c=3.5e-7 and 3.5e-9 (free
  d_max=0.98476 and 0.999922 respectively), but G_c=3.5e-8 stalls on unloading
  at nominal step 21. The failures have finite arrays, nondecreasing accepted
  histories, and source-equation residuals around 1e-7 at final rejection.
  No universal convergence claim is warranted.
- 2026-10-20T19:55:41Z (`date -u`): Manufactured studies extended to n=64.
  Structured c L2 error falls from 6.65e-4 (n=8) to 1.05e-5 (n=64).
  Distorted c falls from 1.47e-3 to 1.93e-4; isolated pressure-gradient error
  remains near 0.025. This limit will be reported explicitly.
- 2026-10-20T20:01:26Z (`date -u`): Began a temporal alternative: lag only
  the history entering the damage equation by one accepted increment, while
  updating returned H/alpha/W from current fields. This is a first-order
  operator split with an exact tangent of its own residual, not an omitted
  tangent block of the fully implicit residual. Host remains monolithic.
- 2026-10-20T20:02:05Z (`date -u`): Lagged alternative completes all three
  formerly stalled physical problems, including the fatigue-floor case
  (H=9.28e19, free d=1). Independent 1D UEL/FV temporal comparisons support
  the same limiting PDE under refinement.
- 2026-10-20T20:04:01Z (`date -u`): Selected the lagged time discretization.
  Archived exact implicit source as verification/alternatives/uel_implicit.f90
  (SHA256 fcd4c37a1a12...). Selected source SHA256 ba0d3006baeb... is identical
  to the lagged comparison source. All cyclic studies are rerun with prefix
  selected_; historical final_* directories denote the earlier implicit
  choice. Updated tests to preload H for equilibrium patches and added an
  analytical three-cycle W/alpha/H test.
- 2026-10-20T20:07:18Z (`date -u`): Selected source passes element tests,
  bitwise rejected-step test, notch and extreme-floor global tests, and all
  n=8/16 D0 checks so far. Normal n=16 D0 takes mean 2.03 Newton updates,
  maximum 3, with no cutbacks. Structured/distorted manufactured n=64 rerun
  with the selected source reproduces the spatial results exactly.
- 2026-10-20T20:11:25Z (`date -u`): Added an independent homogeneous scalar
  implicit-root reference to isolate fatigue time error. Initial reference
  used previous d as a bracket endpoint; roundoff on unloading put the root
  a few ulps outside it. Corrected to the guaranteed physical bracket [0,1].
- 2026-10-20T20:12:54Z (`date -u`): Homogeneous fatigue errors confirm first
  order, but substantial coarse-step bias: max transient d error 0.239 at
  40 steps/cycle, 0.0139 at 640. Alpha errors also decrease approximately by
  two per step halving in the fine range. This is a known accuracy limit.
  Selected large-step D0 completes 20 distorted-mesh cycles directly with
  at most three Newton updates per increment.
- 2026-10-20T20:14:28Z (`date -u`): Controlled selected-time comparison with
  full damage integration still overshoots (d_max=1.15981) versus lumped
  d_max=0.51420. Keep lumped reaction. Fine cyclic runs continue; no further
  formulation changes are planned.
- 2026-10-20T20:19:18Z (`date -u`): Built a separate Q2 verification family,
  with direct displacement Hessians in grad(g*p), consistent original damage
  reaction, and 5x5 quadrature. The generic host is used unchanged with
  reference-only nine-node macroelements; delivered topology remains Q4.
  Reference tangent FD error 6.51e-11; a quadratic-displacement pressure-flux
  patch error is 1.14e-17. Started spatial and temporal reference studies.
- 2026-10-20T20:22:32Z (`date -u`): Q2 n=4,8,16 complete five cycles.
  Cycle-five D_int approaches 0.0644; selected Q1 n=8,16,32 approaches this
  from 0.1636,0.0922,0.0718. This is useful evidence for the pressure/source
  formulation's PDE limit. Q2 n=32 and finer-time checks remain running.
- 2026-10-20T20:27:35Z (`date -u`): Added direct analytic fatigue-floor H
  and nonfinite-input rejection checks; selected element passes both. No
  changes to the selected Fortran source (still ba0d3006baeb...).
- 2026-10-20T20:31:49Z (`date -u`): All eleven selected host CLI studies
  complete, including structured n=64 through five cycles and n=8,16,32
  through twenty cycles on structured/distorted meshes. Structured/distorted
  n=16 time studies at 40,80,160 steps/cycle all pass without cutbacks.
- 2026-10-20T20:33:04Z (`date -u`): Artifact audit passes: every selected
  manifest matches the delivered source, all saved fields/observables are
  finite, H/alpha snapshots are monotone, d stays in bounds, rollback is
  bitwise exact, and all host file hashes match the earlier provenance
  snapshot. Final higher-order comparisons and tables are pending.
- 2026-10-20T20:39:32Z (`date -u`): Finalization complete. All selected cyclic and
  higher-order reference studies finished. Q1 fields and total flux approach
  the Q2/direct-Hessian family on mesh refinement, with reference self-
  convergence recorded. Time-refinement field changes approximately halve
  on both structured and distorted meshes. Final results have no pending
  runs; source and host hashes are unchanged. Updated STATUS.md to record
  completed evidence and remaining accuracy/validation limits.
