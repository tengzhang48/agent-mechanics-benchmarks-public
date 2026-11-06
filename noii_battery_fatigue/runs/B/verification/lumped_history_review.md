# Independent addendum: lumped reaction and history integration

This addendum supersedes any implication in `independent_review.md` that
element consistency checks establish global robustness of the fully implicit
variant. Both consistent and nodally lumped fully implicit damage reactions
failed the distorted n=32, seed=3 notch case despite cutbacks. An active-H
equality tangent steered an earlier n=16 path successfully but did not remove
the later failure. There is no basis for a universal tangent-only cure.

## Reaction extension and its derivative

The inspected extension in
`tmp/lumped_psilag_candidate/element/lump_reaction.py` implements

    R_i = ell^2 sum_j L_ij d_j + m_i d_i - a v_i h_i,
    m_i = sum_G w_G N_i(G),   h_i = sum_G w_G N_i(G) H_phase(G),
    a = 1-kappa.

Here `w_G` includes the physical Jacobian, and the intact constraint remains
the consistent mass equation M[(d-1)+v]=0. Thus v=1-d nodally at an exact
assembled root. Its evaluation order retains a small independently resolved
v when d rounds to 1.

The exact damage-row derivative is

    dR_i/dU_j = ell^2 L_ij [U_j is d]
                 + m_i [U_j is d_i] - a h_i [U_j is v_i]
                 - a v_i sum_G w_G N_i(G) dH_phase(G)/dU_j.

The extension multiplies all four generated H derivative blocks (F, c, d,
vfrac) by the row-node factor -a*v_i and adds precisely the two explicit
nodal diagonals. H has no q or gradient dependence. The retained dummy d
argument ensures that the generated d block exists even when its carrier
derivative is zero. I found no missing product-rule term or sign error.
The assert-guarded text transformation fails rather than silently matching
an unexpected generated layout. The public generator and host are unchanged.

This is a different spatial quadrature from the original consistent reaction,
not an exact rearrangement of that discrete reaction. For smooth fields on
shape-regular meshes it is a consistent nodal product approximation; the
observed O(h^2) damage error supports that assertion. Constants in this
statement need not be uniform as H becomes spatially singular or extremely
large. The small intact field remains necessary to represent the near-one
limit accurately.

For frozen nonnegative H, condensation gives

    [ell^2 L + diag(m_i+a h_i)] d = a h_i.

This matrix is positive definite. A discrete maximum principle additionally
requires the appropriate nonpositive off-diagonal stiffness entries. It is
not guaranteed for all distorted Q4 meshes: independently assembled global
L had 74 positive off-diagonal entries for the host's distorted n=32,
seed=3 mesh (maximum .0847996952); the corresponding structured mesh had
none. Bounds observed in notch tests are useful evidence, not a proof for
every allowed distorted mesh. The same qualification applies to pointwise
nodal irreversibility.

## The full-H lag and its controller

Using accepted H_n in the damage equation while committing physical H_{n+1}
is a first-order splitting on regular, smooth solution paths. It is not the
fully implicit finite-increment equation. It leaves each accepted damage
field out of equilibrium with its just-committed H and can produce an
O(1) relaxation as dt tends to zero from a fixed accepted state. Quasistatic
softening may also have jumps or multiple branches, where a smooth-path
consistency argument alone does not establish convergence to a unique
solution of the stated model.

The concern is observable even in the independent spatially uniform scalar
problem. At 40960 steps/cycle, full-H lag produced 72 sign changes of dW
before t=.8, versus 2 for the implicit scalar reference. Accumulated positive
variation is sensitive to these small oscillations: alpha at that time was
9.39545e-9 versus 9.09446e-9. The final-alpha error was not monotone with
refinement. Merely showing that W itself is close or that H/alpha are
monotone does not control the error in accumulated positive W variation.

An unconditional cutback test based on

    1/(1+a H_old) - 1/(1+a H_trial) <= .05

is not a sound universal remedy. It can reject intermediate Newton trials
before the host line search tries a smaller update. More decisively, a
converged-state test found an actual refinement deadlock. Starting from the
controller's accepted t=.95 notch state, the uncontrolled lumped full-H-lag
element converged with no cuts at all four dt values below:

| dt | maximum history proxy | maximum nodal damage change |
|---:|---:|---:|
| .00625 | .0542182960 | .0592924883 |
| .0001 | .0506767286 | .0592924883 |
| .00001 | .0506606950 | .0592924883 |
| .000001 | .0506593331 | .0592924883 |

The responsible integration point had H_old=2.02248008 and limiting
H_new=2.56894320. Alpha/alpha_T changed from .225451 to .245658, below the
fatigue threshold. Thus the .05 criterion cannot pass this state simply by
time refinement. See `check_lag_controller_deadlock.py` and
`lag_controller_deadlock_review.json`. A fatigue-active gate avoids this
specific event, but would still need evidence that its remaining criterion
does not obstruct a fatigue-driven relaxation.

## Bounded energy-numerator lag with implicit fatigue

The newer candidate instead uses

    psi_phase = min(psi_current, max(psi_old, theta*psi_current)),
    H_phase = max(H_old, 2 ell psi_phase/[Gc f(alpha_current)]),
    alpha_current = alpha_old + max(g(v_current) psi_current-W_old, 0).

It commits the full physical H based on psi_current, plus physical alpha,
W, and psi_current. This retains current damage/fatigue feedback while
lagging the energy numerator that produces the difficult mechanical
softening coupling. All four quantities are per-point SVARS; no persistent
state was found in the wrapper or generated source.

Since 0<=theta<=1, H_phase<=H_physical. At successive accepted steps,

    H_phase[n+1] >= H_physical[n] >= H_phase[n].

Thus the coefficient driving damage is itself irreversible. On smooth paths
with positive energy, psi_old differs from psi_current by O(dt), and the
theta branch becomes inactive as dt tends to zero. This gives a defensible
first-order split. It does not give a convergence proof at damage jumps,
an unconditional Newton guarantee, or the same finite-step result as the
fully implicit model. A nonzero initial energy requires compatible initial
equilibrium or a separately documented initialization procedure.

Theta=.01 is a numerical stabilization parameter, not a material parameter.
It supplies some current-energy feedback on the first/new-energy step; theta=0
allows severe coarse-step fatigue overshoot. It must be declared and tested.
In the scalar study, theta=.001,.01,.1 become indistinguishable by 320
steps/cycle and converge with approximately first-order final-alpha error.
For theta=.01, errors relative to the fine implicit reference are .0952%,
.6088%, .6719%, .3059%, and .00257% at 20,40,80,320,40960 steps/cycle.
These nonmonotone coarse results mean a single favorable coarse error is
not a time-convergence claim. The floor-onset time is much more sensitive
than final alpha and must be reported separately.

`psilag_phase_oracle.py` independently checks the compiled candidate with
different old states at each point. The four points simultaneously select
accepted-energy, theta-current predictor, current-energy cap, and frozen-H
branches, during loading, unloading, and the fatigue floor. It checks
closed-form damage residuals and exact tangent actions for affine strain and
uniform c, vfrac, d, and q variations. Maximum scaled residual, tangent,
and physical-state errors were 1.90e-15, 7.22e-15, and 1.28e-14.
Erasing old psi changed the phase residual by .13269973 while leaving the
physical H/alpha/W/psi output bitwise unchanged, as required. Discarded trial
calls also left the repeated RHS, tangent, SVARS, and PNEWDT bitwise unchanged.
This check used generated source SHA256
`90b72e74a82b697696bbd5ab1ee5f1c6fea2f87d7789bb554dbcb7a2066f88b0`
and wrapper SHA256
`6cb542ed98964c38f7e72a20597c0cc8a7d0a3c4ca815c2e31a5f23df6e6379e`.

## Independent reference scope and recommendation

`uniform_time_reference.py` is independent of the material and generated
assembly. Its scalar bracketing is sound: for fixed prescribed energy,
H_phase(alpha) is nondecreasing, W(alpha) is nonincreasing, and
alpha-alpha_old-max(W-W_old,0) is strictly increasing. The frozen-old-H work
bounds its root. The same reasoning holds for bounded psi lag because its
energy predictor is independent of alpha. Its finest implicit run is a
numerical temporal reference, not an exact continuum solution. The 655360
steps/cycle implicit reference has final alpha=2.8057798063e-7 and floor
onset t=.9223144531; an 80-step onset at 1.0 must not be called exact.

`manufactured_history_oracle.py` is a valid independent spatial test. Its
chosen d=.4+.1 cos(pi x)cos(pi y) satisfies the stated elliptic damage PDE
with its analytically constructed positive frozen H and zero normal flux.
Zero u and c=c0=0 give exactly zero new energy. Five-point error integration
is independent of the element's 2x2 rule. Both consistent and lumped variants
show approximately second-order L2 and first-order gradient convergence on
n=8,16,32,64 structured and distorted meshes. This verifies spatial
consistency of the reaction choice, not the coupled history integrator.

My qualified preference is the resolved intact field with the nodal reaction
extension and bounded energy-numerator lag, retaining implicit current
fatigue. It has a stronger consistency argument and better scalar evidence
than full-H lag with a hard H-proxy controller. Promotion still needs the
new candidate's distorted notch/spatial-floor runs, full coupled time
refinement, and theta sensitivity. No evidence here establishes arbitrary
mesh damage bounds, unique quasistatic evolution, hidden-problem success,
or Abaqus execution.

## Update: bounded phase relaxation

The unrelaxed bounded-psi candidate subsequently passed the distorted n=32
notch, but failed the spatial extreme-fatigue n=16, seed=3 case at step 8,
with c/d residuals about 1.27e-4/2.14e-4 after ten cutbacks. Thus the
preference above is an integration direction, not certification of that
unrelaxed source.

The proposed additional reaction is

    sum_G w_G N_i(G) gamma_old(G) [d_i-d_old_i],
    gamma_old = eta*min(H_old,1000), eta >= 0.

Its coefficient and old nodal damage must be frozen during the increment.
The exact tangent adds only the damage diagonal
sum_G w_G N_i gamma_old. Computing gamma from current H instead would
require additional derivatives and would reintroduce the abrupt floor
coefficient the cap is intended to avoid. U-DU is the host's accepted old
nodal field; a finite-difference check must keep that previous field fixed.

This term can be written tau_n (d_i-d_old_i)/dt with
tau_n=gamma_old*dt<=1000*eta*dt. For fixed eta and cap, tau_n tends to zero
under time refinement. It is vanishing numerical relaxation, not a new
fixed material viscosity. On smooth solution paths the residual perturbation
is O(dt). If damage has bounded total variation, its time-integrated
absolute perturbation is bounded by gamma_max*dt_max*TV(d), which also
vanishes. These facts do not establish the same solution branch at a
discontinuous quasistatic jump: finite per-step relaxation can select a
different limiting jump path even while its integrated residual tends to
zero. The challenge's reference comparison therefore still requires
observable/history refinement through the transition.

For frozen H the condensed matrix gains a nonnegative diagonal and remains
positive definite. On meshes satisfying the required M-matrix conditions,
0<=d<=1 and monotone growth can be propagated from a previous subsolution
with nondecreasing H_phase. Those conditions still need the geometry
qualification above. The extra diagonal does not change the signs of
off-diagonal Q4 stiffness entries and hence does not create a universal
discrete maximum principle. It also cannot by itself certify the coupled
mechanical/concentration Schur complement.

The uniform formula is exactly

    v_new = (1+gamma_old*v_old)/(1+gamma_old+a H_phase).

For fixed H, its error relative to the equilibrium intact fraction contracts
by rho=gamma_old/(1+gamma_old+a H_phase). Because H_phase>=H_old,

    rho <= eta/(eta+a).

For eta=.1,.5,1,2, this upper bound is approximately .091,.333,.500,.667.
The bounded old-H choice therefore gives a finite per-step memory fraction,
and its memory becomes negligible when H_phase is at the fatigue floor
far above 1000. Using d_old rather than an independently rounded 1-v_old
changes the last few bits near d=1; with bounded gamma this is a small
absolute perturbation, unlike multiplying a rounded gap by unbounded H.
The intact constraint and its physical degradation parity still need
explicit reporting.

The scalar alpha root remains monotone. Its numerator and gamma are frozen,
the denominator grows with alpha, and squared v hence decreases with alpha.
The correct upper work bound uses v=(1+gamma*v_old)/(1+gamma+a H_old).
The updated independent scalar reference implements these expressions.
For eta=2, theta=.01, the final-alpha errors at 320,1280,10240 steps/cycle
are .447%, .165%, .0202%; onset times .75,.86875,.91640625 approach the
fine implicit reference .92231445. This supports a vanishing error but
does not justify calling coarse onset times accurate.

A bounded, declared, time-refined relaxation is preferable to a known
spatial extreme-fatigue failure if it completes the relevant coupled
tests. The smallest eta with demonstrated global robustness is a reasonable
choice. It should be accompanied by eta/time sensitivity and the explicit
statement that finite-step damage is relaxed and jump-path agreement is
not proved. The whole-case global result remains separate from the extension
review below.

### Compiled relaxation extension check

The inspected viscosity extension computes gamma immediately after reading
the accepted GP H, uses the generated nodal `d_old=U-DU`, and adds exactly
the nodal reaction and positive damage diagonal above. Its coefficient is
not recomputed from trial H. The wrapper requires 14 props and 16 SVARS,
validates nonnegative eta/cap, and retains local trial-state commit behavior.
No persistent or external state was introduced.

`psilag_phase_oracle.py` now accepts `--element-root` (default the workspace
root), compiles that root without using its candidate helper scripts, and
also covers the viscosity. The previous nodal damage values are deliberately
different at all four nodes. Besides the three mixed-branch history cases,
a fourth case activates the old-H cap while keeping all coefficients
moderate enough to resolve the added reaction in floating point.

The candidate generated source SHA256 was
`bd65169e050795991f890d6b485875e49abaed037b53934681600b6a11fd7b6e`
and wrapper SHA256
`f74b5c1d58279c747f738b54bb83dc0996eacfe58cbf3ba3041e42096be2a931`.
At eta=.5, cap=1000, the full residual/tangent/history checks passed with
the same approximately 1e-14 scaled errors as the unrelaxed source.
Independently setting eta=0 isolated an added matrix that differs from the
exact gamma diagonal by at most 4.27e-14 in absolute value, and an added
reaction differing by at most 8.89e-15. Physical H/alpha/W/psi outputs are
bitwise unchanged by this parameter change.

Changing only the accepted previous nodal damage, leaving U and all old
SVARS fixed, gives the predicted residual change -Gamma_i*delta(d_old_i)
within 3.34e-15. Its tangent and physical state output remain bitwise
identical. This detects accidental use of current nodal damage or an
interpolated old GP value in place of old row-node damage. Mixed psi-history
branches, discarded-call repeatability, and the accepted-psi removal control
continue to pass. Results are in `psilag_visc_phase_oracle.json`; rerun
against the promoted final root to bind this evidence to its delivered hash.

### Whole-reaction and current-phase relaxation

Whole-case tests subsequently rejected the assumption that passing a restart
establishes a robust viscosity parameter. The old-H-only coefficient failed
the spatial extreme-fatigue case at eta=.1 and .5. Whole accepted-reaction
scaling instead uses

    Gamma = eta [1+min(a H_old,Hcap)].

It supplies a positive initial mass scale even when old H is tiny and gives
the cleaner uniform memory bound rho<=eta/(1+eta). It retains bounded
vanishing tau=Gamma*dt and an exactly frozen coefficient. Its independent
compiled phase oracle also passed (source `627c093d...`), but whole-case
failure/tiny-step evidence means that mathematical consistency is not a
global convergence guarantee for this variant either.

The next candidate uses the current **phase** carrier:

    Gamma(H_phase) = eta [1+min(a H_phase,Hcap)].

The exact additional H-chain factor in the damage row is now

    -a v_i + Gamma_H (d_i-d_old_i),
    Gamma_H = eta*a below the cap, zero above the cap.

The explicit damage diagonal remains 1+Gamma and the explicit intact
diagonal remains -a H_phase. These are the terms implemented in the
inspected current-phase extension. All current H derivatives go through the
generated H_phase carrier; none use the separately committed physical H.
That distinction matters: their values differ in the bounded-psi split.

Before the cap, uniform condensation gives exactly

    v_new = eta/(1+eta) v_old + 1/(1+eta) [1/(1+a H_phase)].

Thus the scalar W(alpha) is still nonincreasing, and the alpha root remains
unique. At the uniform root the H-chain coefficient simplifies to
-a/(1+aH), while the damage restoring diagonal is multiplied by 1+eta.
Above the cap the ordinary fixed-Gamma monotonicity proof applies. The
coefficient is continuous at the cap, with a semismooth derivative switch.
Away from a root, however, the H-chain factor can change sign; the added
term is not an unconditional positive contribution to the full nonlinear
Jacobian. Neither uniform uniqueness nor an added diagonal proves global
coupled Newton convergence on distorted meshes.

The current-phase source
`5b1a9102d3ace2ec5894f68e3a1f054d0cc3b4abcd1c7f0d618a3adcb8cbf0f2`
passed the independent compiled oracle. It now forms the complete analytic
added Gamma tangent over all 24 DOFs, including Gamma_H, and isolates that
matrix by subtracting an eta=0 call. Maximum absolute error was 2.85e-14;
the isolated added reaction error was <=1.96e-14. Changing only accepted
old nodal damage predicts both the reaction and Gamma_H tangent changes
within 3.43e-14. Physical H/alpha/W/psi outputs remain bitwise unchanged.
Results are `psilag_currentstiff_phase_oracle.json`.
The oracle accepts an explicit `--relaxation-form` for the three candidates;
this selects an independent mathematical expectation rather than inspecting
the model source to infer its behavior.

If all fixed-eta variants still fail, a distinct principled alternative is
eta_eff=C/sqrt(dt) in the current-phase expression. Its physical viscosity
is bounded by C(1+Hcap)sqrt(dt), so it vanishes while its per-step damping
actually increases on cutback. A fixed accepted state then approaches a
small damage increment as dt tends to zero, rather than retaining a finite
algorithmic relaxation. This predicts O(sqrt(dt)) smooth temporal error,
slower than fixed-Gamma first order, and still needs jump-path and positive
W-variation refinement. Very small dt can also create large cancellation
between the bounded-H reaction and the viscosity reaction; its floating-point
range must be tested rather than assumed. This is a fallback proposal,
not a checked delivered method.
