# Formulation and evidence

The delivered element uses six Q1 nodal fields on four-node quadrilaterals,
with a continuous pressure projection, an intact-fraction variable, a lumped
damage reaction, and a bounded lag of the tensile-energy numerator in the
phase equation. Accepted physical histories use the stated current-state
formulas. This split has a finite-time-step error; it is not the fully
implicit discrete system. The independent uniform reference shows convergence
of that error on fine time grids. A severe heterogeneous fatigue-floor test
still fails, so general nonlinear robustness has not been established.

The physical conventions follow PUBLIC_TASK.md exactly: small strain,
2D identity/trace/deviator, p_eff=tr_2(sigma_eff)/3, and
K=lambda+2mu/3. Thickness and all parameters are dimensionless. No outside
paper, hidden problem or independent full coupled reference was available.
The generator carries 3x3 tensors, but only the active upper 2x2 stress block
enters assembly; psi and pressure use explicit 2D formulas. The unused third
stress component does not introduce a plane-strain constitutive law.

## Spatial discretization

Unknowns per node are ux,uy,c,d,q,vfrac, mapped to Abaqus DOFs
1,2,11,12,13,14. Geometry and all fields are isoparametric bilinear Q1.
The host solves one unsymmetric monolithic system. With a=1-kappa,
v=vfrac, g=kappa+a*v^2, N_i the scalar basis, B_i=grad(N_i), and R=-RHS,

    Ru_i = sum_GP w grad(test_u_i) : g sigma_eff
    Rc_i = sum_GP w [N_i (c-c_old)/dt
                    + B_i.(g D grad(c) - chi c grad(q))]
    Rd_i = sum_GP w [ell^2 B_i.grad(d)
                    + N_i (d_i - a v_i H_phase)]
    Rq_i = sum_GP w N_i (q - g p_eff)
    Rv_i = sum_GP w N_i ((d-1) + v).

c_old is interpolated from U-DU. Concentration storage is backward Euler.
Only the damage reaction uses the row node's d_i,v_i instead of their GP
interpolants. Histories, constitutive laws and other terms use GP values.
There is no damage clipping, added viscosity, penalty on
concentration, or changed natural boundary condition.

### Pressure gradient

q is the global consistent L2 projection of g p_eff into continuous Q1.
Its positive-definite mass block determines all q nodal values, including
the constant mode, without a gauge or boundary condition. Transport uses grad(q), including both
strain-driven and degradation-driven variation of degraded pressure.
Elementwise broken second displacement derivatives would lose jumps in
strain: for example I_h(x^2) has zero broken second x derivative on every
rectangular refinement, although the exact derivative is 2.

Projection is a recovery choice with limitations. An exactly representable
bilinear displacement/linear-pressure patch has the correct nonzero flux and
q-u tangent. Independent constrained-swelling transport has second-order
concentration and first-order flux convergence. A constant-effective-pressure
test with varying damage separately checks p_eff*grad(g).

These examples do not prove strong convergence of every recovered gradient.
The independent recovery test ux=I_h(x^2), q=L2(du_x/dx)
finds an O(h^1/2) gradient boundary layer on structured meshes. On randomly
distorted meshes the potential converges approximately as h but its gradient
has persistent fine-scale oscillations: the unnormalized exact q=2x test
has grad(q) L2 errors about .54 at n=128 for seeds 1/3. See
verification/pressure_projection_refinement.json for the precise convention.
Thus a generic claim of strong distorted-mesh pressure-gradient convergence
would be false. Coupled field/observable convergence must be measured
separately; this is a material remaining limitation of the recovery.

The frozen unconverged severe-floor diagnostic in
verification/pressure_floor_diagnostic/results.json finds large q/g
amplification near the degradation floor even though q nearly satisfies its
L2 projection equation. q/g is a diagnostic quotient only; the UEL contains
no division by g. The alternative product g*L2(p_eff) gives a 14.77-times
larger global drift norm, with the largest gradient discrepancies outside
the lowest-g region. Neither recovered gradient is an exact reference.
The causal role of pressure recovery in failed Newton, and a validated
remedy, remain unverified.

The actual steady transport oracle with prescribed ux=.01*x^2, uy=0 and
d=0 provides a further check. Its exact zero-flux solution satisfies
D*log(c/c0)+chi*B*(c-c0)=chi*P*x, where
P=4*.01*(mu+lambda)/3 and B=4*beta*(mu+lambda)/3.
On n=8,16,32,64, concentration L2 error decreases from .001165 to .000162
on structured meshes and .001103 to .000171 on distorted seed3 meshes;
the last refinement orders are .980 and 1.124. At n64 the C_int error is
about 1.40e-5 on both meshes. Distorted concentration/pressure gradients
still oscillate, while flux L2 error decreases from .000903 to .000254.
This demonstrates concentration and integral convergence for this prescribed
displacement problem despite the strong-gradient limitation; it does not
establish convergence for every coupled fracture path. See
verification/quadratic_drift_results.json.

### Intact fraction and lumped reaction

Rv=M(d+v-1)=0 implies v=1-d nodally in real arithmetic because the consistent
mass matrix is invertible and every v DOF is free, including next to
prescribed damage. The extra field preserves an intact gap smaller than an
ulp of d. For uniform fixed H,

    d=aH/(1+aH),     v=1/(1+aH).

At H=1e15, the nearest representable uniform d leaves a scalar reaction
residual about 8.09e-4 when evaluated using 1-d. The independent v representation reduces the
uniform residual to roundoff. Compiled global solves from the host's
initial d=v=0 pass through H=1e20. A cyclic uniform test actually reaches
f=1e-15 and H=9.35e17. v is corrected by its equation during the first
increment; no undocumented nonzero initial auxiliary values are assumed.

The consistent GP reaction still suffers cancellation for strongly
heterogeneous H inside a cell: interpolated v may combine large nodal
contributions when multiplied by huge H. The selected nodal reaction avoids
that cancellation. For fixed H it is

    ell^2 L d + diag(m) d - a diag(S) v = 0,
    m_i=sum_GP w N_i,       S_i=sum_GP w N_i H.

After eliminating v, the reaction diagonal is positive. This is a different
finite-mesh discretization, consistent for smooth fields. An independent
Neumann elliptic problem with analytic d=.4+.1*cos(pi*x)*cos(pi*y) gives d
L2 errors .002282,.0005792,.0001454,.00003637 for structured n=8,16,32,64;
distorted errors are .002445,.0006624,.0001652,.00004159. Gradient errors
converge approximately as h. A deliberately scaled H source is detected.

The delivered source also passes independent analytic small-matrix checks
with H=1e15 at selected GPs, H=1 at others, and a two-element history jump.
Intact-fraction relative errors are below 1e-15 without cutbacks. These
frozen-coefficient checks isolate representation and reaction conditioning;
they do not remove the evolving heterogeneous-fatigue failure.

The earlier fully implicit parent formulations additionally allow a same-time
comparison of consistent versus lumped reaction. At the first D0 cycle end,
the structured max-nodal damage difference decreases from .007728 at n16
to .001642 at n32; the damage-integral difference decreases from 1.405e-4
to 1.620e-5. This is supplementary spatial-discretization evidence from
historical sources, not final-source temporal validation. The exact records
are in verification/alternatives/implicit_lumped/verification/.

Positive lumped reaction does not guarantee a maximum principle for arbitrary
Q4 geometry: distorted stiffness matrices can have positive off-diagonals.
No general nodal bound or damage-irreversibility proof is claimed. Observed
nodal bounds, decreases and v=1-d error are recorded separately from monotone
GP histories. No clipping hides them.

There is a useful restricted comparison argument on a fixed, axis-aligned
square structured mesh. Assume a>=0, H_phase>=0, positive quadrature weights,
damage boundary values fixed in [0,1], and zero natural damage flux on the
remaining boundary. The Q1 basis is nonnegative at
the Gauss points, so m_i>0 and S_i>=0. L has nonpositive off-diagonals and
zero row sums. Thus A=ell^2*L+diag(m+a*S) is a strictly diagonally dominant
M-matrix with nonnegative inverse. With natural damage conditions, A*d=a*S
and A*(1-d)=m give 0<=d<=1; the same bounds follow after Dirichlet
elimination because its free principal matrix is also an M-matrix.
Between successive converged roots with the same boundary values, the free
rows satisfy A_new*(d_new-d_old)=a*(S_new-S_old)*(1-d_old). Pointwise
H_phase monotonicity below supplies S_new>=S_old, so this right-hand side
is nonnegative and d_new>=d_old. Here A is evaluated at each root; the
argument does not establish existence, uniqueness or Newton convergence
for the nonlinear coupled problem. It applies in real arithmetic between
equilibrated states and does not extend to arbitrary distorted cells or
remove algebraic roundoff.

## Histories and temporal split

At each GP, H0,alpha0,W0,psi0 are committed start-of-increment values and
remain fixed throughout Newton and line search. The physical trial update is

    psi = (K/2)*max(tr(ee),0)^2
          + mu*((ee_xx-ee_yy)^2/2 + 2*ee_xy^2)
    W = g(v)*psi
    alpha = alpha0 + max(W-W0,0)
    f = 1                                      if alpha<=alpha_T
        max(1e-15,max(0,1-k_fg*log10(alpha/alpha_T))^2) otherwise
    H_phys = max(H0,2*ell*psi/(G_c*f)).

The phase equation alone uses

    psi_phase = min(psi, max(psi0, theta*psi)),    theta=.01
    H_phase = max(H0,2*ell*psi_phase/(G_c*f(alpha))).

Fatigue remains implicit in current alpha and v. On unloading the upper cap
uses the current psi; on increasing psi only its numerator is lagged. The
small predictor floor supplies a current-energy drive on the first nonzero loading step
when psi0=0. theta is declared property12. Committed output is always
H_phys,alpha,W,psi, never H_phase. H_phase<=H_phys; at the next accepted
increment H_phase>=the previous physical H. For smooth positive psi and
small dt the predictor floor is inactive and numerator lag is O(dt).
This supports first-order consistency on regular paths, not a convergence
proof across fracture jumps or general coupled bifurcations.

Independent scalar quasistatic bracketing solves the uniform physical model
without importing the UEL or generator. At theta=.01 its final-alpha
relative errors versus the fine fully implicit scalar reference are .306%,
.0809%,.01026%,.002569% at 320,1280,10240,40960 steps/cycle. These fine grids
recover the reference's two work-direction reversals before t=.8.
Coarse errors are nonmonotone: at 80
steps/cycle the final-alpha error is .672% and floor onset is t=.6 versus
about .9223 in the fine reference. Successful Newton convergence is not
an accuracy guarantee. Actual uniform FEM results agree with the independent
split scalar calculation at their sampled accepted times.

A separate 20-cycle D0 temporal study reaches sustained fatigue. Runs
through 1280 steps/cycle complete without cutbacks, but have not established
asymptotic temporal convergence. The exact piecewise-linear D_int relative
time L2 differences are 1.927% for 320 versus 640 and 1.849% for 640 versus
1280; corresponding maximum curve differences are 9.80% and 6.53%.
The interpolated first d_max=.9 times on 40/80/160/320/640/1280 grids are
159.53/143.74/135.70/135.09/134.70/127.55. Final damage L2 differences
between successive 80->160->320->640->1280 grids increase from .001403 to
.002033/.002421/.002735. Reflection asymmetry at cycle16 is below 1.5e-8,
so it does not explain the material path difference there.
One additional 2560-step run is pending within the organizer budget; no
further formulation change accompanies it. These results limit time-accuracy
claims from successful solves. verification/long_time_refinement/ retains
the source-matched interim records, including the historical trapezoidal
320-grid comparison and the exact-polyline 1280-grid comparison.

There are 16 SVARS in cyclic GP order (--,+-,++,-+) at +/-1/sqrt(3).
Each block is H,alpha,W,psi. H slots are 1,5,9,13 and alpha slots 2,6,10,14.
All H/alpha slots are declared in the contract. W and psi may decrease.
Only the accepted converged-state call is retained by the host. The wrapper
uses a temporary state array and copies it out after finite valid evaluation
on LFLAGS(3)=1. There is no SAVE, COMMON, element cache, file I/O, or mutable
module state. A forced rejection before any UEL evaluation versus after a
converged trial commit produces bitwise identical accepted results.

## Alternatives and remaining failure

- The initial five-field fully implicit formulation loses absolute residual
  accuracy near complete damage. Log-H and gap-scaling alternatives were
  tested; a current-gap scaling can falsely converge at d=1 when the correct
  d is about .5. These are rejected, with artifacts retained.
- Adding v fixes uniform high-H roundoff. Fully implicit H nevertheless
  fails a restrained, precracked side-charging case at an old-H switch.
  Exact derivatives, alternate equality choices, and a smaller complex-step
  size do not remove a nondescending local Newton direction. A tie policy
  that helps n16 still fails n32. Lumping alone does not repair this.
- Lagging all H makes several difficult cases solve, but an independent
  scalar test exposes spurious work reversals on fine grids and
  positive-variation bias persisting through 655360 steps/cycle. It is not
  the delivered method.
- Hard H/alpha increment caps can deadlock under cutback because the
  quasistatic phase relaxation remains finite as dt shrinks. Several bounded
  algorithmic relaxation variants also fail full paths despite exact local
  tangents. Their sources and outcomes are labeled alternatives.

The selected split improves the tested notched path and retains fine-time
scalar convergence. A harsher authored case combines spatially varying
history with alpha_T=1e-10; it fails at nominal step8 after cutbacks, even
though the uniform test reaches the same fatigue floor. The element is
therefore not demonstrated robust for all heterogeneous near-complete-damage
paths. STATUS.md gives final completed source-matched studies and failures.

## Checks, quadrature and generator provenance

Full 2x2 Gauss quadrature integrates the mass constraints exactly but not all
nonlinear terms or distorted gradient rational functions. A 3x3 regeneration
through the public ElementConfig API measures sensitivity. Field residual
characteristic scales are .01 for ux/uy,1 for c/d and .001 for q/v; they
change convergence weighting, not equations.

All twenty final-source primary public D0 cases complete with zero cutbacks:
five-cycle structured/distorted n=8,16,32,64 mesh runs, fixed-mesh time
refinement through 1280 steps/cycle, and three twenty-cycle runs. The 145
cycle snapshots are finite and their independently recomputed host
observables agree within 4.81e-12 relative error. This bridge audit verifies
the saved outputs, not their accuracy against an independent PDE solution.

At cycle5, damage L2 differences to structured n64 are .03950 and .008194
for structured n16/n32, and .04743/.009997/.0005342 for distorted
n16/n32/n64. The maximum D_int curve difference between the two n64 meshes
is .628% of the structured reference peak. Structured n32 versus n64 still
differs by 9.35% in that metric, so the finest mesh is a reference for
self-convergence, not a proven mesh-independent answer. The fixed-n16
five-cycle 640-versus-1280-step D_int curve difference is .469%.

Matched twenty-cycle structured runs at eighty steps/cycle give final
(t=160) d_max=0.991334757457 and D_int=0.147650168259 at n16, versus
d_max=0.365011371816 and D_int=0.081797971433 at n32. Both runs use the
same frozen source and time grid. The large difference shows residual
spatial sensitivity at this fixed time step; it does not establish a
mesh-independent late-cycle solution. Joint fine space/time convergence
through the late damage transition remains unestablished. The six-grid
sustained-fatigue time study isolates dt at fixed n16 and is reported
separately above.

The standard-signature wrapper initializes work arrays and checks dimensions,
finite inputs, material ranges and corner Jacobian positivity. Invalid
states or DTIME<=0 return finite zeros, unchanged SVARS and PNEWDT=.5.
Unsupported procedure requests return initialized zeros without committing.
The stock det(F)>0 trial guard is retained conservatively. The supported
contract is quasistatic small strain, one RHS, no distributed loads.
The physical model has no explicit time dependence apart from storage and
accepted history. The generator's TIME(2)=0 initialization branch is not
entered during the host's positive-end-time solves; DTIME<=0 is rejected
before the generated core is called.
No Abaqus execution is claimed.

The bundled abaqus_ufl at organizer-declared commit a02ed7c translates authored
Python physics and supplies interpolation, quadrature, physical-history
updates, complex-step material/H-carrier tangents and most assembly.
No host or package/template file was changed. The generator does not supply
this fatigue model or the nodal reaction scheme. lump_reaction.py makes a
checked deterministic output extension: it replaces the H-carrier phase
reaction, multiplies its generated derivatives by -a*v_i and adds exact d/v
nodal diagonals. Both entry names are changed for the checked wrapper.
Independent assembled residual and FD checks test the extended result;
the stock material verify gate alone would not do so.

The final Fortran has 4055 unchanged generated lines, 22 rewritten/added core
lines and 90 handwritten wrapper lines: 4167 total, 97.31% unchanged generated
by this conservative line-count rule. Generic unused templates inflate the
percentage; it does not measure authorship of physics. generation.json and
provenance.json record counts, hashes, compiler identity and unchanged
host/tool audits. No independent full coupled reference or hidden-problem
execution is available; mesh/time/quadrature comparisons are self-convergence
evidence and must be read with the pressure-recovery limitation above.
