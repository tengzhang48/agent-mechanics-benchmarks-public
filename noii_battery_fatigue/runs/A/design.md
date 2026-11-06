# Direct UEL for chemo-mechanics and fatigue damage

This is an offline arm-A implementation written from PUBLIC_TASK.md. It uses
the supplied host unchanged, exclusively for compilation/calling and global
solving. No Abaqus execution, external code, generators, or external literature
was used. Numerical evidence and limitations are in `verification/` and
`STATUS.md`.

## Fields and equations

There are six isoparametric bilinear fields per node, in this order:
`ux, uy, c, d, q, z`. The two additional fields are

\[
q=g(d)p_{\rm eff},\qquad z=A(1-d),\qquad A=(1-\kappa)H.
\]

In the selected first-order time split, A in the damage/source equations is
evaluated from the accepted start-of-increment H_n. Current H is still updated
and returned from the current trial fields. The equation above denotes the
limiting continuous model; the temporal discretization is explained below.

The additional fields require no prescribed boundary values. Neither has a spatial derivative
in its own constitutive equation. A derivative of q appears in the physical
concentration flux. For positive kappa and finite H, eliminating q,z in the
continuous equations recovers the specified model exactly.

All strains use the two-dimensional symmetric displacement gradient. Swelling
subtracts `beta*(c-c0)` from each of exx, eyy. The Lamé constants have the stated
three-dimensional formulas, while the stress trace has two components and
pressure divides by **three**. The history split is exactly

\[
\psi^+=\tfrac12(\lambda+2\mu/3)\langle e_{xx}+e_{yy}\rangle_+^2
 +\tfrac\mu2(e_{xx}-e_{yy})^2+2\mu e_{xy}^2.
\]

In particular, deviatoric compression contributes. The mechanical stress uses
the full elastic energy, with the stated g degradation; it is not replaced by
a tensile-only stress. Unit thickness is used.

## Spatial discretization and stress-driven diffusion

All six fields use Q1 interpolation on the actual nonaffine map. The Jacobian
is recomputed at each integration point; physical gradients use its inverse
transpose. A tensor-product 3x3 Gauss rule integrates the constitutive and
history terms. This also integrates the degree-four products in affine Q1
mechanics/degradation more accurately than 2x2 integration. On nonaffine maps
the integrands are rational; this rule is not asserted to be exact.

With test function N_a, the unsymmetric residuals are

\[
 R_{u_a}=\int B_a^T\sigma,\quad
 R_{c_a}=\int N_a(c-c_n)/\Delta t+
       \nabla N_a\cdot[gD\nabla c-\chi c\nabla q],
\]
\[
 R_{q_a}=\int N_a(q-gp_{\rm eff}).
\]

The q equation is a **global consistent L2 projection**, solved together with
the other fields. It is not an elementwise interpolation of four local stress
values. Its assembled mass block is positive definite. There is no pressure
boundary condition and no artificial zero pressure-gradient condition. The
concentration equation supplies the physical zero total flux boundary
condition. Both the g-gradient and effective-pressure-gradient contributions
are contained in grad(q). Backward Euler uses the supplied `U-DU` for c_n.

Why not differentiate the Q1 stress inside each element? For example,
u_x=x²/2, u_y=0, c=c0, d=0 on rectangular meshes gives a true pressure
proportional to x. The interpolated displacement has piecewise constant
dux/dx; its broken pressure gradient is identically zero at every refinement.
`verification/test_reference.py` demonstrates this inconsistency and the
convergence of the projected alternative. Ignoring this part of the flux or
retaining only grad(g) would solve a different problem.

Why not simply smooth q with a fixed Helmholtz length? That introduces another
physical length and boundary conditions. A vanishing smoothing length could
be a valid alternative, but would need its own accuracy study. It was not
implemented here. A mixed volumetric-stress alternative, using
`sigma_xx=g*mu*(exx-eyy)+1.5*q` and the corresponding yy stress, was tested.
It slightly improves some pressure errors but does not remove the distorted
mesh gradient noise. The selected method retains direct constitutive stress.

**Observed limit:** the smooth manufactured problem gives approximately
second-order L2 errors for u,c,d,q and first-order flux/grad(q) on structured
meshes. On randomly distorted meshes u,d still approach order two, while c,q
approach order one and flux decreases more slowly. The isolated grad(q) error
stays around 0.024--0.026 through n=64. Thus this formulation does **not** have
demonstrated strong H1 convergence of pressure on those meshes. The total
flux and primary fields converge in that test; cancellation of oscillatory
concentration/pressure gradients matters. This is a material limitation, not
a claim of optimal recovery on arbitrary quadrilaterals.

A separate higher-order 2D comparison provides a further check: affine Q2
displacement/concentration/damage, direct displacement Hessians in grad(g*p),
the original consistent damage equation, and 5x5 quadrature. It has no q,z
auxiliary fields. It uses reference-only nine-node macroelements on the
unchanged generic host; the delivered element remains four-node. Reference
macro n=4,8,16,32 and time refinement are checked, including its own tangent
and quadratic-displacement flux patch.

Against Q2 n=32 at cycle five, selected structured Q1 d L2 differences are
0.123,0.0350,0.00954,0.00281 on n=8,16,32,64; c L2 differences decrease
0.0184,0.00528,0.00135,0.000358. Flux L2 difference falls from 0.00880 to
0.000354 (the reference flux norm is 0.0330). Distorted Q1 fields and total
flux also approach that reference, while isolated grad(q) remains noisy.
The Q2 n=16 versus n=32 d difference is 0.000206, substantially smaller than
the Q1 differences. This supports convergence of the chosen formulation;
it does not prove accuracy of every hidden crack/localization problem.

## Damage equation and saturation

Directly forming `(1+A)*d-A` loses precision for large A. The fatigue floor can
make H exceed the reciprocal of machine precision, so even writing
`d-A*(1-d)` cannot represent a small complement at d=1. Introduce z and use

\[
 -\ell^2\Delta d+d-z=0,\qquad
 i z+s(d-1)=0,\quad i=(1+A)^{-1},\ s=A/(1+A).
\]

The second equation is the original algebraic source relation divided by a
strictly positive coefficient. Its residual and derivatives stay bounded as
H grows: the limiting equations impose d=1 and determine finite z from the
first equation. This is tested at H=1e22. No H cap, fatigue modification or
damage clipping is used.

The initially tried fully integrated Q1 d,z source split overshot d=1 on
coarse D0 meshes (n=8 maximum 1.154 over 20 cycles). The selected formulation
uses nodal quadrature of the **d,z factors** in their zero-order terms:

\[
 R_{d_a}=\ell^2\int\nabla N_a\cdot\nabla d+m_a(d_a-z_a),
\quad R_{z_a}=\int N_a[i z_a+s(d_a-1)],\quad m_a=\int N_a.
\]

H, i and s are still evaluated at nine Gauss points;
history is not averaged or transferred between increments. After global
assembly let b_a=integral(N_a*i). Eliminating z gives a mass-lumped reaction
with effective coefficient `A_a=m_a/b_a-1`. This is a positive, harmonic
average of `1+A`. It is consistent as the mesh refines and nondecreasing if
every stored H is nondecreasing. On rectangular meshes the diffusion
stiffness has the nonpositive off-diagonal entries needed for a discrete
comparison principle. That sign property is not guaranteed on every skew
quad, so no universal nodal-bound theorem is claimed for arbitrary geometry.

The manufactured d error remains order two, with a larger constant than the
fully integrated alternative. The independent transient FV study verifies
that the harmonic-average/lumped method approaches the original damage PDE.
The coarse implicit-history/lumped D0 maximum was 0.5283. Separate distorted
implicit tests reached free damage 0.999922 without overshoot. The selected
lagged-history version also reaches full saturation in the tests below.
With the selected time split, the fully integrated alternative reaches
d_max=1.15981 on n=8 over 20 cycles, while the selected lumped version reaches
0.514202. Thus the benefit of lumping survives the temporal pivot.
Small accepted nodal damage decreases can remain at global-solve accuracy;
H and alpha never decrease.

## History, time integration and tangent

Each of nine points stores a triplet `(H,alpha,W)`, ordered with xi varying
fastest, eta slowest. SVARS are 27 doubles; `history_slots` lists all H and
alpha indices. The initial values are zero. Each call computes, from the
incoming start-of-increment SVARS,

1. `W_trial=g(d)*psi_plus`;
2. `alpha_trial=alpha_n+max(W_trial-W_n,0)`;
3. the exact clipped fatigue factor and `G_f=G_c*f`;
4. `H_trial=max(H_n,2*ell*psi_plus/G_f)`.

The returned W is the **current** W, not its running maximum. This is what
allows alpha to accumulate over repeated unloading/reloading. Only the
host's accepted final call commits these current history values. A prescribed
three-cycle strain test agrees with analytical alpha, H and W at every point.

The selected damage balance uses **H_n**, not H_trial. This is a first-order
operator split in time: solve the mechanics/concentration/damage equations
with the accepted history, then advance the history from the accepted new
fields. There is one host Newton over all six fields, with a block-triangular
history dependence, and no outer staggered iterations. At smooth times the
lag is O(dt); it vanishes under time refinement. This is a time discretization
error, not an exact solution of the endpoint algebraic damage equation at
finite dt. It can delay damage growth and affect fatigue at coarse steps.

The fully implicit candidate used H_trial and differentiated all active
history branches. It failed some notch/localization cases despite cutbacks.
The selected split completes those cases, including a saturated-fatigue
case with H about 9e19. Independent finite-volume comparisons use both the
same time split (to isolate spatial error) and an implicit-history reference
under dt refinement (to check the limiting PDE). These comparisons, rather
than an assertion of accuracy at large dt, justify the split. Exact implicit
source is archived as `verification/alternatives/uel_implicit.f90`.

No SAVE variables, files, COMMON blocks, caches, element-id tables or clock
heuristics exist in the UEL. TIME/KSTEP/KINC are not used to guess whether a
call is accepted. The host owns every transaction. Rejected calls and line
searches only produce disposable trial SVARS. A global forced rejection test
before versus after high-history trial calls is bitwise identical in U,
SVARS and observables.

The analytic tangent includes every coupling in the selected time-discrete
residual. H_n, i and s are fixed during Newton, so the source derivative is
`i*dz+s*dd`, with nodal d,z factors for the selected quadrature. Active alpha,
fatigue and max branches affect the returned history, not the current damage
residual. The concentration tangent differentiates both c and grad(q), g,
and the backward-Euler mass term. It is intentionally unsymmetric.

The returned sign is `RHS=-R`, `AMATRX=dR/dU`, including off-diagonal blocks.
Finite-difference tests keep U-DU fixed and test loading, truly frozen H and
alpha, and the fatigue floor. The smallest relative matrix errors are around
1e-11. Analytic nonaffine patches also check the compression/deviatoric law,
equibiaxial tension, rigid motion and mechanical force integration.

Backward Euler and endpoint positive-increase accumulation need time-step
refinement: a step can miss an interior W peak. Both D0 and the independent
transient problem are studied at successively halved dt. No higher-order
temporal claim is made.

A homogeneous fatigue benchmark isolates this limitation without spatial
error. An independent scalar bracketed root solves the fully implicit
history/damage law at 5120 steps/cycle, with a 2560-step self-check. With
G_c=3.5e-6 and alpha_T=1e-6, the selected scheme's maximum transient d error
is 0.239, 0.117, 0.0572, 0.0280, 0.0139 at 40,80,160,320,640 steps/cycle.
The corresponding alpha error also decreases approximately at order one.
This demonstrates consistency, but also shows that a converged Newton solve
at 40 steps/cycle can have substantial early fatigue error. Extreme-floor
survival tests establish finite behavior, not resolved temporal accuracy.

## Host contract and failure behavior

The contract has six nodal DOFs, Abaqus numbers `[1,2,11,12,13,14]`, all eleven
named model parameters, nine history points, and `[1,0,1,0,0,0]` LFLAGS. The
characteristic residual scales are 0.1 for ux,uy,c,q and 1 for d,z. These are
fixed and do not suppress failed residuals. `JPROPS(1)=1` selects lumping;
zero reproduces the fully integrated comparison. An optional JPROPS(2)=1
selects the mixed volumetric stress experiment; absent/zero uses the final
direct stress law.

The UEL initializes finite zero work arrays. Unsupported dimensions, invalid
material values, nonfinite inputs, extremely large trial values, or a
nonpositive/nearly singular integration-point Jacobian request `PNEWDT=0.5`
and leave SVARS untouched. DTIME<=0 has the same finite rejection behavior.
The element computes into local arrays and only copies SVARS after all
points and finite-output checks pass. The normal dt>0 path has no arbitrary
step cap; Newton failures trigger the host's transactional cutbacks. ENERGY(2)
is the integrated degraded full elastic energy; other energy entries are
zero and no fracture dissipation accounting is claimed.

This is not a proof of unconditional nonlinear robustness. The old implicit
candidate completed large-step D0 through cutbacks, but failed some stronger
notch/localization and intermediate-toughness cases. Nonfinite arrays were
not the cause. The selected split completes those documented cases, including
free d=1, while its coarse-time accuracy requires refinement.
The host has not been modified to stagger fields, change its line search,
increase iteration limits, or relax tolerances. Hidden problems and actual
Abaqus operation remain unverified.
