# Independent formulation review

For the delivered bounded-energy-lag source, see
[final_independent_review.md](final_independent_review.md). This file retains
the earlier fully implicit snapshots and their explicitly stated hashes.

Review date: 2026-10-20 UTC. This review read `ORGANIZER_NOTE.md`,
`PUBLIC_TASK.md`, `host/README.md`, and the bundled abaqus-ufl skill/API/source.
It attacks the proposed mathematics independently; it is not a proof of the
hidden problems or a claim of Abaqus execution.

Final review snapshot: the promoted implementation has six nodal fields
`ux,uy,c,d,q,vfrac`, using the auxiliary intact fraction described below.
Its generated source SHA256 is
`744b5303512e26f261cd5d3158f85e2d8492fb20712bc2a901cb66062812e1ca`;
the checked wrapper SHA256 is
`5608c8db5306f92725c8192fd4f23dfc9d0ab8f459577d4de37d48b48a05dcd9`.
The original d-only high-H subtraction problem identified in this review was
addressed by that field. Final independent compiled results appear at the end.
An additional notch variant later exposed a nonsmooth global Newton failure
of this snapshot. Its localized diagnosis and the subsequent fully implicit
equality-policy candidate are documented in `notch_failure_review.md`.

## Assessment

The proposed extra continuous Q1 nodal field q is a defensible choice:

\[
 (N_a,q-g(d)p_{\rm eff})=0,\qquad
 R_c=(N_a,(c-c_n)/\Delta t)+(\nabla N_a,gD\nabla c-\chi c\nabla q).
\]

It defines a global L2 projection, rather than taking an elementwise second
derivative of Q1 displacement or differentiating an arbitrarily selected side
of an interelement stress jump. The q-q block is the positive-definite mass
matrix. Its constant is determined by the projection, so q requires no
Dirichlet constraint or gauge. This is not an incompressibility constraint and
does not need the u-p incompressibility inf-sup argument.

With valid geometry, kappa > 0, and fixed histories, mechanics retains an
ellipticity lower bound proportional to kappa, scalar diffusion has coefficient
g D >= kappa D, and damage has a coercive reaction-diffusion operator
ell^2 L + M_{1+(1-kappa)H}. These facts do not prove global existence,
uniqueness, positivity of concentration, or convergence of the full nonlinear
fatigue system. The full Jacobian is unsymmetric.

## Highest-priority issues

1. **Projection gradients need nonuniform and mesh-refinement evidence.**
   H1 convergence of q is stronger than L2 convergence of q. Even the smooth
   tensor-product example u_x=x^2, u_y=0, c,d constant, on square grids gives a
   boundary layer in the projection of the broken Q1 derivative. Therefore
   claiming first-order global grad(q) convergence without measurement would
   be incorrect. The independent calculation below gives only order 1/2.
   At moderate distortion the Jacobian and physical gradients vary within an
   element; do not replace them with affine geometry or center-point values.

2. **History trial branches require consistent derivatives.** Old H, alpha,
   and W are fixed arguments throughout every call of an increment. Compute
   trial alpha = alpha_n + max(W_trial-W_n,0), then fatigue f and trial H.
   Differentiate W_trial and these selected trial branches. In particular,
   freezing H in the damage tangent while evaluating trial H in its residual
   is not an exact tangent. Compare centered finite differences only away
   from a switching surface; use one-sided derivatives at a branch boundary.

3. **Monotone histories do not establish nodal damage irreversibility.** In
   the continuum, a comparison principle with Neumann data and increasing H
   gives monotone damage in [0,1]. A full-integration Q1 reaction-diffusion
   matrix need not be an M-matrix, especially on distorted/coarse meshes or
   with a dominant reaction coefficient. Thus nodal overshoot or slight
   damage decreases are possible without decreasing any stored H. Do not
   silently clamp d to [0,1] or to its previous nodal values: that changes the
   specified nonlinear equations unless implemented as a declared constraint
   formulation. Quantify any actual overshoot.

4. **The fatigue floor creates a floating-point and conditioning limit.**
   f = 1e-15 can create H around 1e15 or larger. In the uniform limit the
   exact damage is B/(1+B), B=(1-kappa)H. Its double-precision representation
   cannot generally make d-B(1-d) smaller than O(eps B). For H=1e15 the
   independent calculation gives residual 8.0926984458e-4; for H=1e20 d rounds
   to 1 and the direct reaction is 1. This is not a divergent physical field:
   the error in d is near machine precision. Nevertheless an absolute Newton
   residual tolerance can fail. Test this explicitly before claiming robust
   operation throughout the fatigue floor branch. Multiplying each element's
   damage residual by a different H-dependent factor changes the assembled
   weak equation; it is not automatically harmless row scaling.

5. **Field residual scaling must be checked at fine h and high d.** q has an
   O(h^2) mass block while its physical magnitude shrinks with g near complete
   damage. A small absolute q residual can hide a large relative error in q
   or its gradient. Consistent block scaling and actual flux errors matter
   more than the raw aggregate Newton residual.

## Exact model conventions and generator traps

The fixed task uses a two-dimensional strain and identity, a two-dimensional
trace, p_eff = tr_2(sigma_eff)/3, and dev_2, while retaining
K = lambda + 2 mu/3. In components,

\[
 \psi^+=\frac K2\langle e_{xx}+e_{yy}\rangle_+^2
  +\mu\left[\frac{(e_{xx}-e_{yy})^2}{2}+2e_{xy}^2\right].
\]

The generator's `core.tensor.dev` is explicitly a 3D deviator and `trace`
includes the zz component, so they cannot be applied indiscriminately here.
The actual model inspected after generation uses a 3x3 thermal strain/stress
**carrier**: only its upper 2x2 block enters the generated 2D momentum loops,
and both psi and pressure use explicit 2D component formulas. Its unused zz
carrier entry need not vanish and is not a physical stress component. This
implementation therefore does not silently introduce 3D thermal physics.
Using its full 3D trace or deviator in a later constitutive rewrite would
change the model and should be rejected by the independent tests.

The generated tuple weak residual is storage*N - flux.grad(N). Therefore
`species_flux` returns **-J_c**, and `phase_flux` returns **-ell^2 grad(d)**.
The UEL returns RHS = -R and AMATRX = dR/dU, exactly as the host expects.
`WeakForm.verify()` checks material methods, not the assembled UEL. A green
default/reference-state check cannot exercise a nonzero pressure gradient.

UEL material states are returned through `stress_PK1`; the generator reads
old SVARS at each integration point and writes trial states. Other material
methods may consume those same old-state arguments and must recompute the
same trial H/alpha. Their residuals must not depend on whichever method was
called first. Never use SAVE/module/file state for accepted-history memory.
The stock generator's finite-deformation det(F)>0 guard is conservative for
this small-strain model and should be described as such if retained.

With Q1 fields and distorted bilinear geometry, 2x2 is exact for the q mass
matrix. It is not exact for every nonlinear coupling: N*g(d)*p_eff is already
degree four per coordinate in an affine element in general. Positive 2x2
weights preserve the useful mechanical/mass coercivity properties, but a
comparison with 3x3 or finer quadrature around sharp damage would bound the
quadrature error. This is an accuracy caveat, not a proof of instability.

## Independent calculations actually run

Command, using the organizer's Python environment and workspace TMPDIR:

```
python verification/independent_math_checks.py
```

The script imports neither the generator nor the element. For p=u_x,x=2x and
u_x interpolated from x^2, it analytically assembles the 1D consistent mass
projection of the elementwise derivative. This is the corresponding 2D
Cartesian tensor-product special case with material prefactors divided out.

| n | L2(q-2x) | L2(grad(q)-2) |
|---:|---:|---:|
| 8 | 1.9387086076e-2 | 3.7976784446e-1 |
| 16 | 6.8545521460e-3 | 2.6864247751e-1 |
| 32 | 2.4234501539e-3 | 1.8995892141e-1 |
| 64 | 8.5681901885e-4 | 1.3432124148e-1 |
| 128 | 3.0293126924e-4 | 9.4979460706e-2 |

The rates approach 3/2 for q and 1/2 for grad(q). The boundary value is
approximately h/sqrt(3), while the exact p(0)=0. The independently assembled
projection still converges; this limits an accuracy claim, not its legitimacy.

A separate scalar hand oracle uses E=10, nu=.28, e_xx=.03, e_yy=-.01,
e_xy=.012, d=.4, alpha_n=W_n=.001, alpha_T=.0005, G_c=.002, ell=.12,
H_n=.3. Its expected values are:

| quantity | value |
|---|---:|
| mu | 3.90625 |
| lambda | 4.97159090909091 |
| K | 7.575757575757576 |
| psi_plus | 0.005765151515151515 |
| p_eff | 0.11837121212121211 |
| g | 0.3600064 |
| W_trial | 0.0020754914424242423 |
| alpha_trial | 0.0020754914424242423 |
| f | 0.47737670427768275 |
| H_trial | 1.4492080900029878 |

These independent numerical values were subsequently matched by the final
compiled element in `uniform_history_oracle.py`, as recorded below.

## Closed-form global drift oracle

`verification/constrained_oracle.py` compiles the delivered contract element
once and uses the unmodified host solver. It clamps u=0 at every node, keeps
damage free, and prescribes c_L=c0+.1 and c_R=c0+.8 only on left/right. Since
c>=c0, strain is isotropic compression and psi_plus=H=alpha=W=d=0.

Writing B=4 beta(mu+lambda)/3=0.4142992424242425 gives

\[
 q=-B(c-c_0),\quad J=(D+\chi Bc)\nabla c,\quad
 F(c)=Dc+\frac{\chi B}{2}c^2,
\]

with exact stationary solution F(c(x))=(1-x)F(c_L)+xF(c_R). The stable positive
root is c=2F/[D+sqrt(D^2+2 chi B F)]. The flux is the constant
(F(c_R)-F(c_L),0).

On Cartesian Q1 meshes, the exact nodal values satisfy the discrete steady
equations: integrating (D+chi B c_h)c_h' exactly gives the endpoint difference
of F. On distorted meshes the numerical solution should converge to this
profile. All q values remain solved unknowns with natural projection
conditions. The script checks the exact-profile residual with DU=0, then
solves from a non-exact linear c profile using dt=1e12, making storage
negligible without modifying host code. Independent quadrature computes
field/gradient/flux errors and boundary flux reactions. Omitting drift would
instead give linear concentration and fails the structured exact-nodal test;
the script also records this wrong-oracle separation.

This script passed with the generated source and standard wrapper available
at review time. Results are in `verification/constrained_oracle_results/`.

| mesh | n | L2(c-c_exact) | L2(J-J_exact) | max nodal c error |
|---|---:|---:|---:|---:|
| structured | 8 | 3.4676446026e-4 | 1.8801771285e-3 | 1.30e-14 |
| structured | 16 | 8.6777268871e-5 | 9.4044026508e-4 | 1.33e-14 |
| structured | 32 | 2.1699734544e-5 | 4.7026426558e-4 | 1.37e-14 |
| distorted, seed 1 | 8 | 3.8403682398e-4 | 1.9570996947e-3 | 1.61e-4 |
| distorted, seed 1 | 16 | 9.8327373592e-5 | 9.9038740085e-4 | 5.55e-5 |
| distorted, seed 1 | 32 | 2.5016525718e-5 | 4.9995317941e-4 | 1.93e-5 |

All six final-source solves accepted one increment without cutbacks, produced
max abs(d)<=7.54e-32, and satisfied the nodal q constraint to <=2.22e-16.
Structured exact-profile
free concentration residual norms ranged from 3.50e-16 to 1.59e-15. Both
integrated boundary fluxes matched the exact 0.15600189393939395 to about
1e-14. Concentration field errors approach second order and flux errors first
order on both mesh families. The concentration nodal exactness on structured
grids does not imply that its bilinear interpolant is the exact continuum
field; the reported field error remains nonzero.

The deliberate broken control kept the original analytical oracle but set
chi=0 only in the compiled element's properties. It was rejected by the same
structured test: nodal concentration error 0.02828028000511973 and exact-state
free concentration residual 0.0037618809859421704. Thus this test actually
detects a removed drift term, not merely a code-versus-itself discrepancy.

This oracle exercises swelling-to-pressure coupling, the unusual 2D trace/3
factor, physical gradients on distorted cells, the drift sign, the projected
q relation, global assembly, and the monolithic solver. Since u is clamped,
it does not certify the displacement-gradient part of pressure or fatigue.

## Additional exact/adversarial tests recommended

- **Uniform fixed H:** u=0, c=c0, q=0, alpha=W=0, old H arbitrary and uniform,
  d=B/(1+B). Require zero residual to the achievable floating-point scale and
  Kdd=ell^2 L+(1+B)M. Sweep H=0,1,1e4,1e8,1e12,1e15. This isolates the large-H
  damage/conditioning problem without any history branch ambiguity.
- **Pressure from displacement:** on affine grids take u_x=a*x*y, u_y=0,
  uniform c,d. Then div(u)=a*y, q=g*2(mu+lambda)*a*y/3 and its nonzero constant
  gradient are represented exactly. On distorted cells this displacement is
  not exactly representable; expect refinement rather than exact patch
  equality. Also test rigid rotation and translation, which leave pressure
  and psi unchanged.
- **Branches:** isotropic tension, isotropic compression, pure shear,
  alpha<alpha_T, alpha>alpha_T, fatigue-floor entry, W_trial<W_n, unloading
  with H_trial=H_n, and new maximum H. Centered tangents must remain within
  one branch. A rejected large trial followed by the same accepted state
  must agree bitwise with a run that never saw the rejected trial.
- **Conservation:** on any valid distorted element, sum all concentration
  residual rows and check it equals integral(c-c_n)/dt, because sum grad(N)=0.
  Sum q residual rows and check integral(q-gp_eff). Use nonconstant c,d,q so
  flux signs cannot remain dormant. On a closed unconstrained-concentration
  mesh, C_int must remain constant to the solved residual tolerance.
- **Near complete damage:** test homogeneous d=1 and d=1-1e-8, then spatially
  varying d with the same mean; verify no 1/g or 1/(1-d) singularities were
  introduced. kappa D=1.4e-6 and kappa E=1e-4 remain nonzero.
- **Invalid/zero step:** make DTIME=0 return finite arrays or a cutback,
  leaving accepted history unchanged; deliberately invert geometry and
  verify finite initialized arrays plus PNEWDT<1. Test old histories with
  nonzero values to expose accidental initialization/reset branches.

Remaining uncertainty: the full coupled fatigue problem still needs actual
host mesh/time studies, near-complete-damage runs, and accepted/rejected
increment checks. This independent review does not replace those studies.

## Follow-up: exact remedies for the high-H residual floor

After inspecting the generated core, checked wrapper, and independently
assembled/finite-difference tests, scalar experiments compared five damage
representations. Every scalar residual used a fixed H, a=1-kappa, A=1+aH,
and exact solution d*=aH/A. A line search accepted only decreasing residuals,
with at most 25 Newton iterations and tolerance 1e-8. These are independent
mathematical experiments, not compiled mesh-level tests of a revised UEL.
Reproduce them with `python verification/high_history_representation.py`;
all scalar traces are saved in `verification/high_history_representation.json`.

**Shared nodal scaling is algebraically legitimate.** If every element
contribution to global row i is multiplied by the same positive function of
shared nodal values, then assembly gives R'_i=s_i R_i, preserving roots.
For a state-dependent scale, the exact derivative is
s_i K_ij + R_i (ds_i/dU_j). This is different from element-specific history
scaling, which weights different summands of the same global row differently.

**A scale based only on the current intact fraction is unsafe.** Taking
s=sqrt((1-d)^2+1e-24) gives approximately F_s=z-Az^2 for z=1-d>0.
Its derivative with respect to d is 2Az-1. The small-gap branch therefore
attracts Newton toward d=1 even when the correct root is far away. At H=1,
initial d=.999, this method accepted d=1.0000000000080642 with residual
8.13e-12, whereas the true d*=.4999974999875. This is a deliberately
adversarial counterexample to accepting this scale alone. Starting d=0 at
H=1e15 also failed the 25-iteration limit because the method initially halves
the intact fraction: residual .888 after 25 iterations.

**A scale frozen from accepted nodal d is a safe small change under a limited
validation envelope.** s=max(1e-12,min(1,abs(1-d_n))) is positive, common to
all contributions of a shared nodal row, and constant during the increment.
Its Jacobian is exactly s K. For a spatially uniform accepted solution,
s approximately equals 1/A_n; monotone H makes s A_trial >= approximately 1,
so scaled residual still controls error in d. The scalar jump H_n=1e6 to
H_trial=1e15 converged in one iteration, with residual 8.09e-10. However a
jump directly from d_n=0 to H=1e15 retains the original 8.09e-4 residual floor.
Cutbacks would need to resolve at least one intermediate accepted state. The
uniform result does not prove this claim on a spatially varying damage front.

**An auxiliary log-history scale solves the scalar floor but adds nonlinear
coupling.** A global Q1 projection s of log(1+aH), followed by the positive
shared-nodal multiplier exp(-s_i), locks the scale to actual H rather than
letting an arbitrary d near 1 suppress its own residual. The scalar equations
exp(-s)(A d-aH)=0 and s-log(A)=0 converged in two iterations even for
H=1e20 from zero initial d,s. For a large jump, the first exact Newton update
gives delta(s)=log(A), delta(d) approximately 1-log(A), so its intermediate
damage can be strongly negative. A full coupled line-search test is needed.

**An auxiliary intact-fraction field is the cleanest exact representation.**
Add continuous Q1 v and the algebraic equation

\[
 R_v=(N_a,(d-1)+v)=0.
\]

Use g=kappa+a v^2, W=g psi, and the damage reaction d-a v H. Keep the
positive gradient term ell^2 grad(d).grad(N), so the host continues to solve
and observe its physical d field. At an exact assembled root,
M(d+v-1)=0. The global consistent Q1 mass matrix M is invertible, hence
v=1-d nodally and pointwise within each element. This holds also at nodes
where the host prescribes d, because v itself has no prescribed boundary
values. The mixed system is exactly equivalent to the original discrete
equations in real arithmetic; v adds precision to the small intact fraction.

For frozen H the damage/intact-fraction tangent blocks are

\[
 \begin{pmatrix}\ell^2 L+M&-M_{aH}\\M&M\end{pmatrix}.
\]

Eliminating delta(v)=-delta(d) recovers ell^2 L+M+M_{aH}. With active
fatigue, g_v=2av, alpha_v=I_loading g_v psi, and
H_v=-H (f_alpha/f) alpha_v on the active-history branch. Thus the reaction
derivatives are r_d=1 and r_v=-a(H+v H_v). Condensation gives
ell^2 L+M_{1+aH+a v H_v}, identical to the original d tangent because
H_d=-H_v. Mechanical, q, concentration, and history couplings transform by
the same chain rule; all uses of g and W must use v consistently.

In floating point, v remains representable even when d rounds to 1. The
original d-only reaction subtracts a large coefficient times an inaccurately
represented 1-d. The mixed reaction instead multiplies that coefficient by
the separately accurate v. Evaluate its constraint as (d-1)+v, retaining
the resolved small v rather than rounding d+v to 1 first.

| fixed H | initial d,v | iterations | final scalar residual norm |
|---:|---|---:|---:|
| 1e12 | 0,0 | 1 | 3.21e-17 |
| 1e15 | 0,0 | 1 | 8.09e-19 |
| 1e20 | 0,0 | 1 | 1.00e-20 |
| 1e15 | 0,1 | 2 | 8.09e-19 |

For initial d=v=0, all reported d values matched the nearest-double exact
d*. The v=0 initial state is inconsistent with its mass constraint but is
exactly the host's default for an additional nodal field; the scalar test
therefore exercises that default. At H=1e20, final d=1 and
v=1.000010000100001e-20. There is no representable d-only gap of that size.

The remaining risks for this option are global sparse conditioning at large
H, Newton behavior when trial fatigue crosses its floor, and controlling the
constraint error. A residual scale such as v=.001 makes that equation visible
to the host's absolute residual test. Always report the actual nodal
norm/max of v-(1-d), and its effect on reconstructed g(d), rather than treating
the extra field as a visualization quantity. Full generated, independent
assembled, rejected-history, and global host checks must be repeated after
adding the sixth field.

## Final six-field compiled re-review

The final Python model and generated core use vfrac consistently in g, W,
stress, q projection, species diffusion, and the damage reaction. The intact
constraint is generated as `(d-1)+vfrac` and has zero gradient flux. Damage
retains ell^2 grad(d).grad(test). The standard wrapper checks 24 element DOFs,
uses a local trial copy of SVARS, and commits that copy only after all finite
and monotone-history checks pass. No SAVE, COMMON, module state, or file I/O
was found in either delivered Fortran source. Only the upper 2D stress block
is assembled. This inspection found no extra physical term introduced by the
auxiliary field.

The final contract scales ux and uy residuals by .01, q by .001, and vfrac
by .001, in the host's characteristic-size convention. This makes the two
auxiliary mass equations visible in the host norm. These are fixed field
scales, not element-dependent changes of the weak equation.

The following independent compiled checks passed against the source hashes
at the top of this review:

1. `python verification/constrained_oracle.py`: all six structured/distorted
   n=8,16,32 cases passed again, including an initial vfrac=0 as supplied by
   the host. The final max nodal abs((d-1)+vfrac) was <=4.44e-16 and the q
   constitutive relation error <=2.22e-16. The corresponding degradation
   parity error max abs(g(vfrac)-g(1-d)) was <=8.89e-16.
   Concentration and flux refinement
   results remain the values tabulated above to the printed precision. The
   deliberately removed-drift control still failed its independent exact
   profile with nodal c error .02828028000511984.
2. `python verification/pressure_gradient_patch.py`: on a 1.3 by .8 rectangle,
   u_x=.03*x*y, u_y=0, c=c0, d=.27, vfrac=.73 give the exact q gradient
   (0,.09462085777698863). The species flux is
   (0,-.008279325055486505), with exact nodal c residuals
   [.005381561286066228,.005381561286066228,
   -.005381561286066228,-.005381561286066228]. The compiled c residuals match
   those numbers bit-for-bit; max q residual is 2.52e-17 and the independent
   closed-form q-u tangent block differs by <=1.67e-16. This exercises the
   displacement-gradient contribution that the constrained-compression oracle
   leaves dormant. Removing displacement while retaining q gives q residual
   .013120758945075757 and is correctly detected. Its intact constraint error
   is exactly zero.
3. `python verification/uniform_history_oracle.py`: all four integration
   points on a distorted element match the frozen hand calculation above.
   Loading, unloading, reloading, and fatigue-floor histories are checked.
   Reloading gives alpha=.003632110024242424 and H=2.1337751328249275;
   the floor gives H=6.918181818181818e14. The largest state error divided
   by max(1,abs(expected)) is 1.28e-14. Erasing old W in a deliberately broken
   reload produces alpha=.004150982884848477, separated from the correct
   answer by .0005188728606060605. Thus old-W accumulation affects the
   quantitative oracle rather than merely its qualitative trend.
4. The uniform-history script makes several arbitrary discarded high/low
   strain calls from the same old SVARS, then repeats its original call.
   RHS, tangent, SVARS, and PNEWDT are bitwise identical. This agrees with the
   absence of persistent state found in source inspection. It checks direct
   element history isolation, not the organizer's complete forced-rejection
   global scenario.

These independent checks certify the stated equations and limits; they do
not certify hidden boundary/load cases or replace mesh/time convergence of
the full cyclic fatigue problem. Global high-H sparse conditioning, spatial
damage-front behavior, long cycles, and forced global rejection still need
the separate host studies. No Abaqus execution was performed.

Subsequent global failures and the reaction/history integration alternatives
are reviewed in [lumped_history_review.md](lumped_history_review.md). That
addendum qualifies the historic six-field results above and reports the
independent compiled bounded-energy-lag oracle and the converged-state
history-controller deadlock test.
