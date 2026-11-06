# Status

Implemented the direct arm-A UEL and contract, with all work confined to this
workspace. The host is unchanged. No internet, installations, external
examples, generators, or Abaqus execution were used.

The selected source uses six Q1 nodal fields, continuous projected pressure,
a bounded auxiliary damage-source equation with nodal reaction quadrature,
backward-Euler diffusion, and a **first-order lag of accepted H in the damage
balance**. Returned H, alpha and W are computed from current accepted fields.
See [design.md](design.md) for the equations and [verification/RESULTS.md](verification/RESULTS.md)
for generated evidence. Source SHA256 begins `ba0d3006baeb`.

Verified:

- Compiled residual/tangent finite differences, analytic nonaffine strain and
  force patches, tension/deviatoric energy, rigid motion, conservation, exact
  three-cycle history, finite zero-time rejection and H=1e22 saturation.
- A global rejected increment gives bitwise identical U, SVARS and observables
  whether rejected before or after the UEL sees high-history trial states.
- All eleven selected host CLI studies complete: n=8/16/32 through 20 cycles
  on structured and distorted meshes; structured n=64 through five cycles;
  and 40/80/160-step time studies on both n=16 mesh types. H/alpha are monotone
  and no candidate run needs a cutback. At n=16,
  mean Newton updates are 2.03, maximum 3, with no cutbacks. The 10-step/cycle
  distorted run also completes 20 cycles, with at most three updates.
- Constrained-notch and extreme fatigue-floor tests complete five distorted
  cycles; free damage reaches one and H reaches 9.28e19, with finite arrays.
- Smooth 2D manufactured fields converge through n=64 on structured and
  randomly distorted meshes. An independent finite-volume coupled reference
  confirms spatial convergence and the first-order temporal limit.
- A scalar implicit reference at 5120 steps/cycle verifies history/damage
  time convergence without spatial error.
- A separate Q2/direct-Hessian 2D formulation converges through macro n=32,
  with its own tangent and stress-gradient patch checks. Selected structured
  Q1 damage L2 differences decrease 0.123, 0.0350, 0.00954, 0.00281 on
  n=8,16,32,64. Displacement, concentration, flux and host observables also
  approach the higher-order result. Q2 self-convergence is recorded. This is
  a self-authored numerical comparison, not the organizer's hidden reference.
- The final artifact audit verifies source hashes, unchanged host files,
  finite saved fields/observables, monotone histories and bitwise rollback.

Known limitations:

- The history lag can have large coarse-step fatigue error. In the homogeneous
  stress test, max transient damage error is 0.239 at 40 steps/cycle, falling
  to 0.0139 at 640. Newton convergence is not an accuracy certificate.
- On random distortions, isolated projected-pressure gradient error does not
  decrease through n=64. Primary fields and total flux do converge in the
  manufactured test, but more slowly for concentration/pressure.
- The old fully implicit candidate fails some localization cases; those
  failed attempts and exact source are retained. The selected split completes
  them. Full reaction integration also overshoots d=1 and is not selected.
- Temporal accuracy of the extreme fatigue-floor survival test is unverified.
  Accuracy of strongly localized/notched fracture propagation against an
  independent reference is unverified. The Q2 comparison covers smooth D0.
  Hidden organizer problems, arbitrary quad families and actual Abaqus
  execution are unverified.

All launched studies are finished. Results, reference self-checks and failed
alternatives are under `verification/`; timings and decisions are in
`research_log.md`. No Abaqus execution is claimed.
