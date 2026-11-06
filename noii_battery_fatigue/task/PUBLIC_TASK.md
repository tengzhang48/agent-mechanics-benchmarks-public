# Robust UEL challenge v2: chemo-mechanics with fatigue damage

**Goal.** Build an Abaqus/Standard user element (UEL) that solves the model below **robustly in a global
analysis**. The model is fixed. **The formulation is yours to choose.** Explore the problem, decide how
to discretize it, and show that your element works on meshes, through damage growth and over many load
cycles.

This is a reproduction and engineering task based on a published model (Noii et al. 2024, lithium
diffusion with swelling and fatigue phase-field fracture). It is not a test of new physics.

## 1. The model (fixed)

The fields are the displacement **u** (2D, unit thickness), the concentration c and the damage d.
Everything is dimensionless and small-strain.

**Constitutive laws**

- ε = sym ∇u,  ε_e = ε − β (c − c0) I
- σ_eff = 2μ ε_e + λ tr(ε_e) I,  μ = E / [2(1+ν)],  λ = Eν / [(1+ν)(1−2ν)]
- g(d) = (1 − κ)(1 − d)² + κ,  σ = g(d) σ_eff
- p_eff = tr(σ_eff) / 3
- J_c = g(d) D ∇c − χ c ∇(g(d) p_eff)

Here I is the 2×2 identity and tr the 2D trace; p_eff divides by 3. These conventions come from the
authors' public script. They are part of the model: you may discuss them, but implement them as stated,
so that results can be compared.

**Balance laws** (in Ω, for t > 0)

- ∇·σ = 0
- ∂c/∂t = ∇·J_c
- −ℓ² Δd + d − (1 − κ)(1 − d) H = 0

Natural boundary conditions are traction-free, zero flux and ∇d·n = 0, except where the host prescribes
values.

**History** (pointwise, irreversible)

- ψ⁺ = (K/2) ⟨tr ε_e⟩₊² + μ dev₂(ε_e) : dev₂(ε_e), with K = λ + (2/3)μ
- W = g(d) ψ⁺
- α(t) = accumulated increase of W, so α grows only when W grows
- f(α) = 1 for α ≤ α_T; otherwise f(α) = max(1e-15, [max(0, 1 − k_fg log10(α/α_T))]²), with
  k_fg = 0.5
- G_f = G_c f(α)
- H(t) = max over s ≤ t of 2ℓ ψ⁺(s) / G_f(α(s))

H, α and W may change only on accepted increments. A Newton iteration or a rejected increment must leave
no trace in them.

**Parameters**

| E | ν | β | c0 | κ | D | χ | ℓ | G_c | α_T |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 0.28 | 0.035 | 0.35 | 1e-5 | 0.14 | 0.25 | 0.12 | *set in the problem file* | *set in the problem file* |

## 2. What is yours to decide

- **How to compute ∇(g p_eff).** Note that p_eff already depends on ∇u. How you discretize this term
  is your decision. Justify it with evidence, including how it behaves under mesh refinement.
- **The discretization.** You choose the interpolation and quadrature for each field, and whatever
  nodal or element-local variables your formulation needs (declare them in the contract).
- **The history update.** You choose how and when H, α and W are updated and stored.
- **Failure behaviour.** You choose how the element behaves on invalid states, and how it asks the host
  for a smaller increment (PNEWDT).
- **Anything else** robustness needs.

The mesh topology is fixed to **4-node quadrilaterals**, which may be non-affine and moderately
distorted. Explain every choice and the alternatives you rejected, with evidence.

## 3. The host

`host/` is the solver part of `abaqus_ufl.fe`, which plays the role of Abaqus/Standard. It calls your
compiled UEL through the standard argument list, with:
- one monolithic Newton per increment, with a line search;
- start-of-increment SVARS at every call, kept only after convergence and restored on rejection;
- cutbacks when PNEWDT < 1 or Newton fails.

**Read `host/README.md` first.** It states every argument the host passes, the time convention, and the
format of `element/uel_contract.json`. In the contract you name your PROPS, so the host can set the
model parameters, and you point to your stored H and α, so the scorer can check them.

`host/run_problem.py` runs the development problem D0 (`host/dev_problem/d0.json`, a square without a
crack under cyclic charging) on structured or distorted meshes, and writes observables computed by the
host from your nodal values. Use it to test your element globally.

## 4. What "robust" means (how you will be judged)

The organizer runs hidden problems on the same host. They are similar to the development problem but
use other boundary conditions, loading, parameters and meshes. No element-level answer key is used,
because formulations differ.
1. **Convergence.** Under mesh and time-step refinement, the fields and observables converge to an
   independent reference solution of the PDE. The observables are the maximum damage, the integrals of
   d and c, and the elastic energy, over the cycles. This is tested on structured and distorted meshes.
2. **History.** H and α never decrease. A forced rejected increment leaves no trace. Results are stable
   under time-step refinement.
3. **Solver behaviour.** Newton converges through damage growth up to d near 1. The element survives
   cutbacks, a large-step run and a direct DTIME = 0 call, and returns no NaN.
   - **Hidden state:** two runs reject the same increment, one after your UEL has seen the trial states
     and one before. They must agree exactly, so all memory must be in SVARS.
4. **Your design.** Whether your formulation choice is justified by evidence.
5. **Honesty.** Whether your claims match what you checked.

## 5. Deliverables

- **`element/`:**
  - the UEL source (standard signature);
  - `uel_contract.json` in the format of `host/README.md` §2: DOF names, named PROPS, the SVARS layout
    with `history_slots`, LFLAGS, and residual scales.
- **`design.md`:** the formulation, the alternatives considered, and why you chose as you did.
- **`verification/`:** your tests and your global studies on the host (mesh and time-step studies), with
  results.
- **`research_log.md`:** times taken from the clock (`date`), failed attempts and pivots.
- **`STATUS.md`:** what works, what does not, and what is unverified.

## 6. Conditions

- **Isolation:** offline, no installs, this workspace only.
- **Environment:** `<env>`.
- **No Abaqus.**
- **Budget:** about 6 hours.
- **No human help.**
- **Arm:** the arm line, with or without `abaqus_ufl`, is given in `ORGANIZER_NOTE.md`.
