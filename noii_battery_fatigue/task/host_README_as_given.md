# The v2 host: how your UEL is called

The host is the solver part of `abaqus_ufl.fe` (a development revision, commit 9ac03e8 of 27 September 2026), vendored in
`hostfe/`. It plays the role of Abaqus/Standard. It calls your compiled UEL through the standard
argument list and solves the global problem. The organizer scores your element on this same host.

## 1. How the host calls your UEL

Every value below is fixed; none is left for you to guess.

| Argument | What the host passes |
|---|---|
| NNODE, MCRD | 4, 2. Nodes are counter-clockwise; COORDS(1:2, a) is node a |
| NDOFEL, NRHS, MLVARX | 4 × `dof_per_node`, 1, NDOFEL |
| U, DU | Current nodal values in your `dof_names` order, node by node. DU = U − (accepted values at the start of the increment). So previous values are U − DU |
| PROPS, NPROPS | Your `props`, with every entry named in `prop_names` that is a model parameter overwritten by the problem's value |
| JPROPS, NJPROP | Your `jprops` (NJPROP = its length; 0 if empty) |
| SVARS, NSVARS | Start-of-increment values at every call of the increment. After a converged increment, the host calls your UEL once more at the converged state, and keeps the SVARS that call returns. A rejected increment restores the previous SVARS. Initial values: your `svars_initial` (default 0) |
| LFLAGS | Your `lflags` (default 1, 0, 1, 0, 0, 0). LFLAGS(3) = 1 asks for both RHS and AMATRX |
| TIME(1), TIME(2) | Both equal the time at the **end** of the increment. This is the abaqus_ufl.fe convention; Abaqus itself passes the beginning |
| DTIME | The increment length. It is always > 0 in solves. The scorer also makes one separate call with DTIME = 0 (§4) |
| KSTEP, KINC, JELEM | 1, 1, element index + 1 |
| PNEWDT | 1 on entry. Return a value < 1 to reject the increment: the host restores SVARS and retries with a step halved, up to 10 times |
| RHS, AMATRX | The host uses R = −RHS and K = AMATRX = dR/dU. K may be unsymmetric |
| ENERGY, V, A, PREDEF, PARAMS, loads | Unused: zero arrays and no distributed loads |

**The global solve.**
- **Newton:** one monolithic Newton iteration over all DOFs per increment, with a line search.
- **Convergence:** judged on the residual, weighted per field by your `residual_scale` (default 1),
  and on the update.
- **Dirichlet conditions** are applied by the host, as described in §3.

## 2. `element/uel_contract.json`

```json
{
  "source": "element/uel.f90",              // or a list of files; compiled in this order
  "compile_flags": [],                       // extra gfortran flags, e.g. ["-ffixed-line-length-none"]
  "dof_per_node": 4,                         // more if your formulation adds nodal variables
  "dof_names": ["ux", "uy", "c", "d"],       // must include these four; name any others yourself
  "abaqus_dofs": [1, 2, 11, 12],             // the Abaqus DOF numbers you would use
  "prop_names": ["E", "nu", "beta", "c0", "kappa", "D", "chi", "ell", "G_c", "alpha_T", "k_fg", "my_flag"],
  "props": [10.0, 0.28, 0.035, 0.35, 1e-5, 0.14, 0.25, 0.12, 1.0, 1.0, 0.5, 0.0],
  "jprops": [],
  "nsvars": 0,                               // your SVARS length per element
  "svars_initial": null,                     // or a list of length nsvars, per element
  "history_slots": {"H": [], "alpha": []},   // 1-based SVARS indices of H and of alpha, one per point where you store them
  "lflags": [1, 0, 1, 0, 0, 0],
  "residual_scale": {"c": 1.0}               // optional characteristic size per named field
}
```

The `//` comments above are explanation only; the real file must be plain JSON.

The model parameter names the host fills in are E, nu, beta, c0, kappa, D, chi, ell, G_c, alpha_T and
k_fg. `history_slots` must point at every stored value of H and of α, because the scorer checks that they
never decrease. The host prescribes values only for ux, uy, c and d; any further nodal variables get
natural boundary conditions.

## 3. Problem D0, the development problem

The problem is in `dev_problem/d0.json` (parameters, loading and boundary conditions). There is no crack.
- **Domain:** Ω = [0,1]².
- **Constraints:** u_x = u_y = 0 at (0, 0); u_y = 0 at (1, 0). Otherwise traction-free.
- **Concentration:** c = c0 + A[1 − cos(2πt/T)]/2 on the top edge y = 1; zero flux elsewhere.
- **Damage:** ∇d·n = 0 on the whole boundary.
- **Initial state:** u = 0, c = c0, d = 0, and zero history.

To run it:

```bash
export PATH=<env>/bin:$PATH
python host/run_problem.py . --mesh structured --n 16            # also: --mesh distorted --seed 1
python host/run_problem.py . --mesh structured --n 32 --steps-per-cycle 80 --cycles 5
```

Each run writes `verification/runs/D0_<mesh>_<n>_s<seed>_spc<steps>/`:
- `observables.csv`: d_max, A_damage, D_int, C_int and E_el per accepted step, computed by the host
  from your nodal values (`observables.py`);
- `fields_cycle_XXX.npz`: nodal values and SVARS at the end of each cycle;
- `run.json`: the solver summary, including `history_min_change`, the smallest change of each
  `history_slots` entry between accepted steps.

Meshes: `meshes.py`, with structured or distorted 4-node quadrilaterals and n a multiple of 8.

## 4. What the scorer adds on top of D0

The scorer uses hidden problems with other boundary conditions, loading and parameters, on structured
and distorted meshes, with the same host and the same contract. One of them drives the damage to 1.
On top of those it runs:
- **History checks:** H and α (your `history_slots`) never decrease.
- **A hidden-state test:** two runs reject the same increment, one after your UEL has seen the trial
  states and one before. They must give identical results, so all memory must be in SVARS.
- **A DTIME = 0 call:** one call at a typical state must return finite arrays or PNEWDT < 1. It must not
  crash.
- **A large-step run:** a quarter of the steps per cycle, completed through cutbacks or directly.
- **Newton statistics** per increment.
