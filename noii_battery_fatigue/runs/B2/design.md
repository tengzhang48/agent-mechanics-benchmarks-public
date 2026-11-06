# Formulation and implementation

The target is `PUBLIC_TASK.md`, including its two-dimensional elastic trace,
pressure divided by **3**, two-dimensional deviator, and positive-energy bulk
coefficient `K=lambda+2*mu/3`. The given constitutive and fatigue laws are
retained. Spatial approximation and time splitting are described below.

## Fields and enhanced strain

A four-node isoparametric quadrilateral carries nodal order
`(ux,uy,c,d,p,s)`, 24 element DOFs. `p` approximates `g*p_eff`; `s` is
independently represented intactness. Nodal displacement, concentration,
damage, pressure and intactness use Q1 interpolation. Integration is 2x2
Gauss, with unit thickness.

Four element-local strain-mode amplitudes are solved and condensed afresh on
every call. With `phi_1=1-xi^2`, `phi_2=1-eta^2`, define

    G_m = grad(phi_m) - <grad(phi_m)>_element,
    E_(i,m) = sym(e_i tensor G_m),
    eps_h = sym(grad u_Q1) + sum_j q_j E_j.

The average uses the same weights and Jacobians as the element. This is an
enhanced, nonconforming strain approximation. The UEL uses `eps_h` for stress,
pressure and history; the host reconstructs its observables from the Q1
nodal displacement. Their finite-mesh difference is measured in the studies.

For fixed current concentration and intactness the local system is linear:

    A_ij = sum_q w_q g_q E_i:C:E_j,
    b_i  = sum_q w_q g_q E_i:C:(eps_Q1-beta*(c-c0)*I),
    q = -A^-1 b.

`enhanced_modes.f90` uses a pivoted four-by-four solve. Positive `g>=kappa`
and positive elastic moduli make the independent local mode system positive.
A relative small-pivot guard handles numerical singularity. The zero-average
gradients preserve affine strain patches on non-affine cells.

The compiled element also passes the exact compatible-linear-swelling patch

    c=c0+a*x,
    ux=beta*a*(x^2-y^2)/2,   uy=beta*a*x*y,

with zero elastic strain on an affine cell. Ordinary Q1 displacement cannot
represent this compatible quadratic field's strain. It produces artificial
energy and fatigue on coarse meshes. The local modes substantially reduce
that bias in D0; the comparison is retained in `verification/exploration/`.

For this exact stress-free field on a uniform square mesh, the host's Q1
reconstruction alone gives nonzero energy
`(beta*a)^2*(3*mu+lambda)*h^2/24` on the unit square, although continuum
energy is zero. `check_energy_bridge.py` verifies this analytical identity.
This observation error still counts in the reported host observable; it is
not subtracted. Reference fields sampled onto the production nodes and
reconstructed with Q1 provide a separate diagnostic of its size.

The matched twenty-cycle n=16 study uses 40 steps/cycle and identical
properties, accepted-history splitting, reaction, guard and tolerances.
Removing only the enhancement hooks makes the generated cores byte-identical.
Ordinary Q1 gives final d_max=.953243; enhancement gives .267728; the finite
independent Q4 reference with the same time split gives .267448. Enhanced
maximum damage agrees to .11%, but damage integral differs by 2.56% and the
host Q1 energy by 50.1%. Peak agreement alone is not full-field validation.
These controls and their provenance are in `verification/enhancement_comparison.json`.

## Pressure gradient and bulk residuals

Pressure satisfies a global consistent L2 projection:

    Rp_a = sum_q w_q N_a [p_q-g_q*p_eff,q].

The mass operator is invertible and introduces neither a pressure gauge nor
an extra pressure boundary condition. Transport uses the Q1 gradient of this
solved pressure:

    J_c = g D grad(c) - chi c grad(p),
    Rc_a = sum_q w_q [N_a*(c-c_old)/dt + grad(N_a).J_c],
    Ru_ai = sum_q w_q grad(N_a)_j * sigma_ij.

`c_old=U_c-DU_c`. Stress and `p_eff` use the enhanced elastic strain. Under
refinement the projection tends to `p=g*p_eff`, including the gradient of
degradation in the stress-drift term.

Direct broken second derivatives of ordinary Q1 displacement were rejected.
On rectangles, the Q1 interpolant of `ux=x^2` has zero elementwise `ux,xx`,
although its exact pressure gradient is nonzero. Ignoring strain jumps loses
that leading term. The independent Qp reference instead uses physical shape
Hessians with order at least two and checks its nonzero stress-gradient flux
against a quadratic analytical example.

Pressure recovery remains a limitation. The raw projected pressure gradient
does not converge strongly in the distorted manufactured family tested. The
complete transport flux and primary fields converge. With the selected
enhancement, concentration converges approximately at first order in that
forced problem, including structured meshes. The enhancement trades some
pressure-recovery accuracy for much better strain/history accuracy in D0.

## Intactness and time discretization

The mass constraint

    Rs_a = sum_q w_q N_a [d_q+s_q-1]

implies nodal `d+s=1`, since the Q1 mass matrix is positive definite.
Degradation is evaluated as `(1-kappa)*s^2+kappa`. Representing `s`
independently avoids cancellation in `1-d` near complete damage. The archived
five-field prototype cannot reduce its residual at H near `1e16`; the final
element resolves intactness below `1e-16` in quantitative compiled/global
tests.

Concentration uses backward Euler with a consistent Q1 mass matrix. This
does not ensure a concentration maximum principle on unresolved meshes. In
a fully damaged beta=0 limit the diffusivity is kappa*D=1.4e-6 and the
diffusion length at t=8 is .00335. Compiled n=8/16 solutions agree with an
independent 1D reduction to 6e-15, but its coarse profiles undershoot c0.
Fine 1D refinement resolves the layer and converges to an analytical
continuous-space backward-Euler reference. This is a measured spatial
resolution limitation, not a claim of concentration positivity. It is a
plausible contributor to the coupled severe-case uncertainty, not a complete
attribution. See `verification/damaged_diffusion_resolution.json`.

Damage uses **accepted start-of-increment
history H_n**. Exact endpoint history is computed from converged current
fields and stored for the next increment. This is a first-order history
split: for smooth evolution, replacing `H_(n+1)` by `H_n` has `O(dt)` error.
Time-refinement evidence is therefore necessary.

Let `L_a=sum_q w_q N_a` and
`M_a=sum_q w_q N_a [1+(1-kappa)*H_n,q]`. The damage row is

    Rd_a = ell^2 sum_q w_q grad(N_a).grad(s) + M_a*s_a - L_a,
    dRd_a/ds_b = ell^2 Kgrad_ab + delta_ab M_a.

This is the negative original damage equation after `d=1-s`, with a row-sum
lumped reaction matrix. Its block is positive definite and linear during
each increment. On uniform square meshes, nonpositive gradient off-diagonals give a
maximum principle, `0<=s<=1`, for damage boundary values in `[0,1]`. With
fixed boundary damage and nondecreasing history, they also give nondecreasing
nodal damage, apart from floating-point error. Distorted cells retain
coercivity; a universal maximum principle is not claimed for arbitrary
distortion or arbitrary rectangle aspect ratios. No clipping is used.

The damage PDE occupies `d` rows and the mass constraint occupies `s` rows.
Prescribed nodal damage therefore imposes the corresponding intactness
through the constraint, without an extra boundary condition.

Fully implicit current-history damage was tested first. It completed D0 but
stalled at history switches in a severe n=32 problem despite ten cutbacks.
Captured fields and global Newton replay diagnostics are retained. Smaller
complex steps, a generalized history tangent, reaction-lumping variants,
nodal residual scaling, larger phase characteristic scales and removing the generic
finite-strain guard did not resolve that failure. Accepted-history splitting
completed the same severe problem. The selected scheme retains the original
strict residual tolerances and quantifies its time-splitting error.

A later controlled test repeated current-H Gauss damage with the four enhanced
strain modes. Its complete condensed tangent passed finite differences to
2.4e-9 normalized error. Nevertheless, the severe n=32 solve failed at step
19 after exhausting cutbacks, with a phase-dominated residual and damage
overshoot to 1.0253. Enhancement alone did not remove the implicit-history
failure. The experimental builder, source snapshots, failed trial and solver
report remain in `verification/exploration/implicit_enhanced*`.

## History and transactions

Each Gauss point stores `(H,alpha,W)`. Point order is counter-clockwise from
`(-1/sqrt(3),-1/sqrt(3))`. H slots are `1,4,7,10`, alpha slots `2,5,8,11`,
and W slots `3,6,9,12`. All initially equal zero.

At each call, always from incoming accepted SVARS, compute

    Wnew = g(s)*psi_plus,
    alpha_new = alpha_old + max(Wnew-Wold,0),
    f = 1                                  if alpha_new<=alpha_T;
        max(1e-15,max(0,1-k_fg*log10(alpha_new/alpha_T))^2) otherwise,
    Hnew = max(Hold,2*ell*psi_plus/(G_c*f)).

The increasing branch deliberately evaluates `alpha_old + (Wnew-Wold)`.
The left-associated expression can decrease alpha by an ulp. A compiled
before/after regression demonstrates the correction.

The generated core returns trial SVARS; the reaction uses a saved copy of
incoming H. The host alone keeps returned state after convergence. Every
residual, line-search and rejected call starts from accepted state. No saved
Fortran state, element-number cache or retained local mode is used. Rejection
before and after actual UEL evaluations gives exactly identical complete
output nodal vectors and SVARS arrays in the rollback test.

## Tangent, interface and failure handling

Generated material tangents use complex steps. Static condensation subtracts
the complete Schur term `K_Xq*A^-1*K_qX`, including concentration, intactness
and pressure-constraint couplings. The damage reaction contributes its exact
diagonal tangent. Central differences verify the complete compiled matrix on
a non-affine cell, including the local solve's dependence on current fields.

The host convention is `R=-RHS`, `AMATRX=dR/dU`, with an unsymmetric matrix.
Characteristic residual scales are `.01` for ux, uy and p, and `1` for c, d
and s. Smaller mechanical/pressure scales retain accuracy near `g=kappa`.
The host's normal monolithic Newton, line search, tolerance and cutback limit
are used throughout.

The wrapper checks dimensions, finite data, parameter ranges and nonnegative
accepted H/alpha. Invalid geometry, a singular local solve, nonfinite outputs
or unsupported requests return finite zero residual/matrix arrays with
`PNEWDT=.5`, preserving incoming SVARS. `DTIME<=0` follows the same rule.
Requests 1, 2 and 5 are supported; only request 1 returns updated trial state.
The generic `det(I+grad u)>0` guard is removed because this small-strain model
uses no deformation-gradient inverse or logarithm. The physical geometry
Jacobian must be positive at the integration points used by the mapping
guard. Arbitrarily warped or concave cells are outside the validated mesh
envelope. ENERGY is unused by the host and remains zero;
reported elastic energy is the host's independent nodal reconstruction.

## Provenance and validation scope

The bundled `abaqus_ufl` package and `host/` are unchanged. The Python model,
standard Quad4 generator, field/state packing, material/history translation,
complex-step engine, shape functions, mapping and bulk assembly are used
directly. The deterministic build renames the entry, removes the finite-strain
guard and inserts the local-strain/Schur calls. `enhanced_modes.f90`,
`accepted_damage.f90` and `uel.f90` are handwritten. The generator does not
natively supply this four-mode condensation or nodal reaction lumping.
Exact line counts, transformations and hashes are in
`verification/provenance.json`.

Evidence includes independent matrix-based element assembly, shear/biaxial/
swelling limits, local-mode finite differences, Fourier diffusion, global
fully damaged affine patches, exact rejection equality, manufactured mesh
refinement and independent high-order Python Qp simulations. The Qp reference
uses neither generated constitutive code nor pressure projection; its Newton
and transactional stepping are provided by the same unmodified generic host.
Hidden-problem/reference results and production Abaqus execution remain
unavailable. See `STATUS.md` and `verification/README.md` for completed runs,
measured convergence and practical limitations.
