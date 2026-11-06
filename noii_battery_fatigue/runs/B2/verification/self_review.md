# Single-agent review

Organizer conditions prohibit subagents and human help. This review was done
by the implementing agent; it is not independent review. Independent means
an algebraically separate implementation/oracle in the tests below.

- Checked the literal public conventions: 2D trace, hydrostatic divisor 3,
  dev_2, K=lambda+2mu/3, positive volumetric energy only, log10 fatigue with
  floor 1e-15, and multiplication by g in both transport terms. Independent
  biaxial/shear cases and intentionally wrong generated controls discriminate
  these choices from plausible incorrect alternatives.
- Checked geometry orientation and J/inverse conventions on a non-affine
  cell. Four local modes have weighted zero-average gradients. Their local
  elasticity uses engineering shear consistently. NumPy weighted least
  squares independently reproduces their residual. Finite differences check
  the full Schur correction, including u, c, s and p couplings.
- Checked 24 DOFs, named properties and the 12 state slots against the host.
  Phase PDE is on d rows and the alias on s rows, so prescribed damage is
  enforced through the mass constraint without an auxiliary boundary value.
  Projection mass matrix has no gauge nullspace.
- Checked that damage uses incoming accepted H, while endpoint H/alpha/W use
  converged current enhanced strain. This is a declared first-order time
  split, not a fully implicit endpoint-history solve. It retains the fixed
  continuum laws only in the time-refined limit. Large-step robustness must
  not be described as large-step accuracy.
- Checked the complete state lifecycle: incoming SVARS are copied before
  generated updates; rejected/invalid calls restore that copy; the host
  controls acceptance. Request 2/5 do not advance history. The forced rejection
  comparison agrees bit-for-bit in all saved outputs and state. No physics
  state persists in SAVE, COMMON, a module, a file or a local-mode cache.
  Generated unused eigenvector helpers contain two DATA-initialized integer
  rotation-index tables; these are read-only constants, not changing state.
- A monotonicity guard caught the actual one-ulp alpha accumulation bug.
  The correction groups the positive energy difference before addition;
  a compiled predecessor reproduces the error, and the same former failing
  global run completes after correction. No unreported history clipping was
  introduced.
- Checked zero-time and invalid-input handling, output finiteness, fatigue
  floor, and near-complete damage. Independent intactness avoids subtractive
  loss of 1-d. Archived five-field control confirms the failure it prevents.
  The det(I+grad u) guard is deliberately absent for this small-strain model;
  geometry det(J)>0 remains required.
- Checked the output bridge: all persisted field indices, connectivity,
  finite states and host observables are audited independently. Pressure
  recovery is not claimed strongly convergent on distorted meshes. Q1 host
  elastic energy and enhanced-strain mechanics energy differ at finite mesh;
  both are explicitly compared, with host energy retained as the observable.
- A separate exact stress-free swelling field validates the host's nonzero
  Q1 reconstruction energy against a closed formula at n=8,16,32,64. It is
  proportional to h^2 despite zero continuum energy. The severe reference
  sampled onto production nodes has a similar reconstruction error; neither
  diagnostic changes or subtracts the reported host energy.
- Checked a fully damaged diffusion limit against an independent 1D reduction
  and an analytical continuous-space backward-Euler reference. Consistent
  mass can undershoot the initial concentration on unresolved low-mobility
  layers. Finer 1D refinement resolves it; only n=8/16 were compiled in 2D.
  This limitation is retained in the status and results, without attributing
  every coupled severe-case error to the same mechanism.
- The physical-space interpolation check initially exposed an overstrict
  inverse-map assertion on wrong candidate distorted cells. Corrected point
  location to accept only converged inverses inside a cell. Known quadratic
  fields now reproduce to about 2e-15 for Q2/Q3/Q4 on both mesh families,
  before any production/reference field comparison is trusted.
- Checked provenance: deterministic regeneration is byte-identical, the four
  final sources are snapshotted before each authoritative run, generator and
  host are not edited, generated-vs-handwritten accounting includes its
  counting convention, and older formulations are filtered by source hash.
  No claim is made that the generator natively supports the manual local
  condensation or reaction lumping.

Finite reference resolution, time-splitting error, pressure recovery,
distorted-cell damage bounds, hidden problems and Abaqus remain explicit
limitations. Numerical tables and final completion status are in results.md
and ../STATUS.md.
