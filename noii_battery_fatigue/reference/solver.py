"""Reference solver for the chemo-mechanical fatigue phase-field model (challenge v2).

Independent answer key (dolfinx 0.10).  Fields: displacement u (2D), concentration c,
damage d and the auxiliary field q = g(d) p_eff, tied to g p_eff by an L2 projection
(int (q - g p_eff) r dx = 0) so that the stress-gradient flux uses grad q.

Model (fixed, see README.md):
  eps_e   = sym grad u - beta (c - c0) I            (I: 2x2 identity)
  sig_eff = 2 mu eps_e + lam tr(eps_e) I,  g(d) = (1-kappa)(1-d)^2 + kappa,  sigma = g sig_eff
  p_eff   = tr(sig_eff)/3
  J_c     = g D grad c - chi c grad(g p_eff)        ->  = g D grad c - chi c grad q
  div sigma = 0;  dc/dt = div J_c;  -ell^2 lap d + d - (1-kappa)(1-d) H = 0
History at quadrature points, updated only after a converged (accepted) step:
  psi+ = K/2 <tr eps_e>_+^2 + mu dev2(eps_e):dev2(eps_e),  K = lam + 2/3 mu
  W = g psi+;  alpha += max(W_new - W_old, 0);  f(alpha) (log fatigue, k_fg);  G_f = Gc f
  H = max(H_old, 2 ell psi+_new / G_f(alpha_new))
The damage equation of a step uses H of the last accepted step (lagged history).
Time: backward Euler.  Nonlinear solver: PETSc SNES Newton (basic line search), LU (MUMPS).
"""
import math
import os
import time
from dataclasses import dataclass, asdict

import numpy as np
from mpi4py import MPI
from petsc4py import PETSc

import basix
import basix.ufl
import ufl
from dolfinx import fem, geometry
from dolfinx import mesh as dmesh
import dolfinx.fem.petsc as fpetsc

HERE = os.path.dirname(os.path.abspath(__file__))
JIT_DIR = os.path.join(HERE, "tmp", "jit")


# ----------------------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------------------
@dataclass
class Params:
    E: float = 10.0
    nu: float = 0.28
    beta: float = 0.035
    c0: float = 0.35
    kappa: float = 1e-5
    D: float = 0.14
    chi: float = 0.025
    ell: float = 0.12
    kfg: float = 0.5
    fmin: float = 1e-15
    # calibrated
    Gc: float = 1e-2
    alphaT: float = 1e-3
    A: float = 0.4
    T: float = 1.0

    @property
    def mu(self):
        return self.E / (2.0 * (1.0 + self.nu))

    @property
    def lam(self):
        return self.E * self.nu / ((1.0 + self.nu) * (1.0 - 2.0 * self.nu))

    @property
    def K(self):
        return self.lam + 2.0 / 3.0 * self.mu

    def c_b(self, t):
        return self.c0 + self.A * (1.0 - math.cos(2.0 * math.pi * t / self.T)) / 2.0

    def asdict(self):
        return asdict(self)


# ----------------------------------------------------------------------------------------
# Constitutive expressions (UFL)
# ----------------------------------------------------------------------------------------
def kinematics(u, c, d, prm):
    I = ufl.Identity(2)
    eps = ufl.sym(ufl.grad(u))
    eps_e = eps - prm.beta * (c - prm.c0) * I
    sig_eff = 2.0 * prm.mu * eps_e + prm.lam * ufl.tr(eps_e) * I
    g = (1.0 - prm.kappa) * (1.0 - d) ** 2 + prm.kappa
    p_eff = ufl.tr(sig_eff) / 3.0
    return eps_e, sig_eff, g, p_eff


def psi_plus(eps_e, prm):
    tr = ufl.tr(eps_e)
    trp = 0.5 * (tr + abs(tr))
    dev = eps_e - 0.5 * tr * ufl.Identity(2)
    return 0.5 * prm.K * trp ** 2 + prm.mu * ufl.inner(dev, dev)


def psi_full(eps_e, prm):
    return prm.mu * ufl.inner(eps_e, eps_e) + 0.5 * prm.lam * ufl.tr(eps_e) ** 2


def fatigue_f(alpha, prm):
    """f(alpha), numpy, pointwise."""
    f = np.ones_like(alpha)
    m = alpha > prm.alphaT
    if np.any(m):
        r = np.maximum(0.0, 1.0 - prm.kfg * np.log10(alpha[m] / prm.alphaT))
        f[m] = np.maximum(prm.fmin, r * r)
    return f


# ----------------------------------------------------------------------------------------
# Geometry helpers
# ----------------------------------------------------------------------------------------
def tri_area_superlevel(v, level):
    """Fraction of a triangle where the linear interpolant of vertex values v (n,3) >= level."""
    s = np.sort(v, axis=1)
    a, b, c = s[:, 0], s[:, 1], s[:, 2]
    frac = np.zeros(len(v))
    frac[a >= level] = 1.0
    m1 = (a < level) & (b >= level)  # level in (a, b]
    with np.errstate(divide="ignore", invalid="ignore"):
        below = (level - a[m1]) ** 2 / ((b[m1] - a[m1]) * (c[m1] - a[m1]))
    frac[m1] = 1.0 - below
    m2 = (b < level) & (c > level)  # level in (b, c)
    with np.errstate(divide="ignore", invalid="ignore"):
        frac[m2] = (c[m2] - level) ** 2 / ((c[m2] - a[m2]) * (c[m2] - b[m2]))
    return frac


def ref_lattice(m):
    """Lattice points of order m on the reference triangle and its m^2 sub-triangles."""
    idx = {}
    pts = []
    for j in range(m + 1):
        for i in range(m + 1 - j):
            idx[(i, j)] = len(pts)
            pts.append((i / m, j / m))
    tris = []
    for j in range(m):
        for i in range(m - j):
            tris.append((idx[(i, j)], idx[(i + 1, j)], idx[(i, j + 1)]))
            if i + j < m - 1:
                tris.append((idx[(i + 1, j)], idx[(i + 1, j + 1)], idx[(i, j + 1)]))
    return np.array(pts, dtype=float), np.array(tris, dtype=np.int64)


class PointEvaluator:
    """Evaluate functions at fixed points (MPI-safe: each point is evaluated on one rank that owns a
    containing cell, values are summed over ranks so every rank gets all values)."""

    def __init__(self, msh, points2d, allow_outside=False):
        comm = msh.comm
        pts = np.zeros((len(points2d), 3))
        pts[:, :2] = points2d
        nloc = msh.topology.index_map(msh.topology.dim).size_local
        tree = geometry.bb_tree(msh, msh.topology.dim, entities=np.arange(nloc, dtype=np.int32))
        cand = geometry.compute_collisions_points(tree, pts)
        coll = geometry.compute_colliding_cells(msh, cand, pts)
        cells = np.full(len(pts), -1, dtype=np.int32)
        for i in range(len(pts)):
            lk = coll.links(i)
            if len(lk) > 0:
                cells[i] = lk[0]
        mine = np.where(cells >= 0, comm.rank, comm.size).astype(np.int32)
        owner = np.empty_like(mine)
        comm.Allreduce(mine, owner, op=MPI.MIN)
        self.outside = owner == comm.size
        if np.any(self.outside) and not allow_outside:
            raise RuntimeError(f"{np.sum(self.outside)} points not found in mesh")
        self.mask = owner == comm.rank
        self.pts = pts[self.mask]
        self.cells = cells[self.mask]
        self.npts = len(pts)
        self.comm = comm

    def __call__(self, f):
        bs = int(np.prod(f.ufl_shape)) if len(f.ufl_shape) else 1
        loc = np.zeros((self.npts, bs))
        if len(self.cells):
            loc[self.mask] = f.eval(self.pts, self.cells).reshape(-1, bs)
        out = np.zeros_like(loc)
        self.comm.Allreduce(loc, out, op=MPI.SUM)
        out[self.outside] = np.nan
        return out


# ----------------------------------------------------------------------------------------
# Solver
# ----------------------------------------------------------------------------------------
def residual_form(w, w_old, tests, H, dt, chi, prm, dx, steady=False):
    u, c, d, q = ufl.split(w)
    _, c_o, _, _ = ufl.split(w_old)
    v, wc, wd, r = tests
    eps_e, sig_eff, g, p_eff = kinematics(u, c, d, prm)
    J = g * prm.D * ufl.grad(c) - chi * c * ufl.grad(q)
    F_u = ufl.inner(g * sig_eff, ufl.sym(ufl.grad(v))) * dx
    if steady:
        F_c = ufl.inner(J, ufl.grad(wc)) * dx
    else:
        # multiplied by dt (row scaling only)
        F_c = ((c - c_o) * wc + dt * ufl.inner(J, ufl.grad(wc))) * dx
    F_d = (prm.ell ** 2 * ufl.inner(ufl.grad(d), ufl.grad(wd)) + d * wd
           - (1.0 - prm.kappa) * (1.0 - d) * H * wd) * dx
    F_q = (q - g * p_eff) * r * dx
    return F_u + F_c + F_d + F_q


class CMFSolver:
    """Monolithic (u, c, d, q) solver with quadrature-point history."""

    def __init__(self, msh, prm, deg_u=2, deg_s=1, qdeg=4, qdeg_obs=None, steady=False,
                 prefix="cmf_", snes_rtol=1e-10, snes_atol=1e-12, snes_max_it=30,
                 linesearch="basic"):
        self.msh = msh
        self.comm = msh.comm
        self.prm = prm
        self.deg_u, self.deg_s, self.qdeg = deg_u, deg_s, qdeg
        self.qdeg_obs = qdeg_obs if qdeg_obs is not None else max(qdeg, 2 * deg_s + 2 * (deg_u - 1) + 2)
        cell = msh.basix_cell()
        Pu = basix.ufl.element("Lagrange", cell, deg_u, shape=(2,))
        Ps = basix.ufl.element("Lagrange", cell, deg_s)
        self.W = fem.functionspace(msh, basix.ufl.mixed_element([Pu, Ps, Ps, Ps]))
        Qe = basix.ufl.quadrature_element(cell, value_shape=(), scheme="default", degree=qdeg)
        self.Q = fem.functionspace(msh, Qe)
        self.w = fem.Function(self.W, name="w")
        self.w_old = fem.Function(self.W, name="w_old")
        self.H = fem.Function(self.Q, name="H")
        self.alpha = fem.Function(self.Q, name="alpha")
        self.Wh = fem.Function(self.Q, name="W")
        self.ind = fem.Function(self.Q, name="ind")
        self.dt = fem.Constant(msh, PETSc.ScalarType(1.0))
        self.chi = fem.Constant(msh, PETSc.ScalarType(prm.chi))
        self.dx = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": qdeg,
                                                          "quadrature_scheme": "default"})
        self.dxo = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": self.qdeg_obs})
        self.tests = ufl.TestFunctions(self.W)
        self.F = residual_form(self.w, self.w_old, self.tests, self.H, self.dt, self.chi, prm,
                               self.dx, steady=steady)
        # collapsed sub-spaces and maps
        self.Vu, self.map_u = self.W.sub(0).collapse()
        self.Vc, self.map_c = self.W.sub(1).collapse()
        self.Vd, self.map_d = self.W.sub(2).collapse()
        self.Vq, self.map_q = self.W.sub(3).collapse()
        self.map_u = np.asarray(self.map_u)
        self.map_c = np.asarray(self.map_c)
        self.map_d = np.asarray(self.map_d)
        self.map_q = np.asarray(self.map_q)
        self.uf = fem.Function(self.Vu)
        self.cf = fem.Function(self.Vc)
        self.df = fem.Function(self.Vd)
        self.qf = fem.Function(self.Vq)
        # history expressions at quadrature points
        u, c, d, q = ufl.split(self.w)
        eps_e, sig_eff, g, p_eff = kinematics(u, c, d, prm)
        X = self.Q.element.interpolation_points
        self.psi_expr = fem.Expression(psi_plus(eps_e, prm), X, jit_options=self.jit())
        self.g_expr = fem.Expression(g, X, jit_options=self.jit())
        self.psi_q = fem.Function(self.Q)
        self.g_q = fem.Function(self.Q)
        # observable forms
        n = ufl.FacetNormal(msh)
        J = g * prm.D * ufl.grad(c) - self.chi * c * ufl.grad(q)
        self.form_D = fem.form(d * self.dxo, jit_options=self.jit())
        self.form_C = fem.form(c * self.dxo, jit_options=self.jit())
        self.form_E = fem.form(g * psi_full(eps_e, prm) * self.dxo, jit_options=self.jit())
        self.form_flux = fem.form(ufl.inner(J, n) * ufl.ds(domain=msh, metadata={"quadrature_degree": self.qdeg_obs}),
                                  jit_options=self.jit())
        self.form_ind = fem.form(self.ind * self.dx, jit_options=self.jit())
        self.snes_opts = {
            "snes_type": "newtonls",
            "snes_linesearch_type": linesearch,
            "snes_rtol": snes_rtol,
            "snes_atol": snes_atol,
            "snes_stol": 0.0,
            "snes_max_it": snes_max_it,
            "ksp_type": "preonly",
            "pc_type": "lu",
            "pc_factor_mat_solver_type": "mumps",
            "mat_mumps_icntl_14": 200,
        }
        self.prefix = prefix
        self.bcs = []
        self.bc_updaters = []
        self.problem = None
        self.newton_log = []
        self.notch_dofs_d = np.zeros(0, dtype=np.int64)  # collapsed-d indices on the notch
        self.cell_cache = {}

    @staticmethod
    def jit():
        return {"cache_dir": JIT_DIR}

    # ---------------------------------------------------------------- set-up
    def finalize(self, extra_F=None):
        F = self.F if extra_F is None else self.F + extra_F
        self.problem = fpetsc.NonlinearProblem(F, self.w, petsc_options_prefix=self.prefix,
                                               bcs=self.bcs, petsc_options=self.snes_opts,
                                               jit_options=self.jit())
        snes = self.problem.solver
        self._res = []
        snes.setMonitor(lambda s, it, rn: self._res.append(float(rn)))

    def apply_bcs(self, t):
        for f in self.bc_updaters:
            f(t)
        for bc in self.bcs:
            bc.set(self.w.x.array)

    # ---------------------------------------------------------------- time stepping
    def _solve(self, t, dt):
        self.dt.value = dt
        self.w.x.array[:] = self.w_old.x.array
        self.apply_bcs(t)
        self._res = []
        self.problem.solve()
        snes = self.problem.solver
        reason = snes.getConvergedReason()
        its = snes.getIterationNumber()
        ok = (reason > 0) and np.all(np.isfinite(self.w.x.array))
        return ok, reason, its, list(self._res)

    def accept(self):
        """History update from converged fields; then w_old <- w."""
        prm = self.prm
        self.psi_q.interpolate(self.psi_expr)
        self.g_q.interpolate(self.g_expr)
        psi = self.psi_q.x.array
        W_new = self.g_q.x.array * psi
        dW = np.maximum(W_new - self.Wh.x.array, 0.0)
        alpha_new = self.alpha.x.array + dW
        f = fatigue_f(alpha_new, prm)
        H_cand = 2.0 * prm.ell * psi / (prm.Gc * f)
        H_new = np.maximum(self.H.x.array, H_cand)
        self.hist_dmin = (self.comm.allreduce(float(np.min(H_new - self.H.x.array)), op=MPI.MIN),
                          self.comm.allreduce(float(np.min(alpha_new - self.alpha.x.array)), op=MPI.MIN))
        self.H.x.array[:] = H_new
        self.alpha.x.array[:] = alpha_new
        self.Wh.x.array[:] = W_new
        self.w_old.x.array[:] = self.w.x.array

    def advance(self, t_old, t_new, depth=0, max_depth=6):
        """Advance from t_old to t_new; on failure bisect (history untouched until acceptance)."""
        dt = t_new - t_old
        ok, reason, its, res = self._solve(t_new, dt)
        if ok:
            self.accept()
            self.newton_log.append(dict(t=t_new, dt=dt, its=its, reason=reason, res=res, depth=depth))
            return [(t_new, its, res)]
        self.newton_log.append(dict(t=t_new, dt=dt, its=its, reason=reason, res=res, depth=depth,
                                    failed=True))
        if depth >= max_depth:
            raise RuntimeError(f"Newton failed at t={t_new} after {depth} bisections")
        self.w.x.array[:] = self.w_old.x.array
        tm = 0.5 * (t_old + t_new)
        out = self.advance(t_old, tm, depth + 1, max_depth)
        out += self.advance(tm, t_new, depth + 1, max_depth)
        return out

    # ---------------------------------------------------------------- field access
    def sync_collapsed(self):
        a = self.w.x.array
        self.uf.x.array[:] = a[self.map_u]
        self.cf.x.array[:] = a[self.map_c]
        self.df.x.array[:] = a[self.map_d]
        self.qf.x.array[:] = a[self.map_q]

    def d_values(self):
        return self.w.x.array[self.map_d]

    # ---------------------------------------------------------------- observables
    def _lattice_data(self, V, m):
        key = (id(V), m)
        if key in self.cell_cache:
            return self.cell_cache[key]
        pts, tris = ref_lattice(m)
        be = V.element.basix_element
        tab = be.tabulate(0, pts)[0, :, :, 0]  # (npts, ndofs)
        gx = self.msh.geometry.x[:, :2]
        gdm = self.msh.geometry.dofmap
        X0 = gx[gdm[:, 0]]
        X1 = gx[gdm[:, 1]]
        X2 = gx[gdm[:, 2]]
        area = 0.5 * np.abs((X1[:, 0] - X0[:, 0]) * (X2[:, 1] - X0[:, 1])
                            - (X2[:, 0] - X0[:, 0]) * (X1[:, 1] - X0[:, 1]))
        # physical x of lattice points: X0 + xi*(X1-X0) + eta*(X2-X0)
        latx = X0[:, None, 0] + pts[None, :, 0] * (X1[:, None, 0] - X0[:, None, 0]) \
            + pts[None, :, 1] * (X2[:, None, 0] - X0[:, None, 0])
        nloc = self.msh.topology.index_map(2).size_local
        dm = V.dofmap.list[:nloc]
        data = (tab, tris, area[:nloc] / len(tris), latx[:nloc], dm)
        self.cell_cache[key] = data
        return data

    def superlevel(self, vals, V, level=0.5, m=None):
        """Area of {f >= level} and max x reached, f = FE function with dof values vals."""
        if m is None:
            m = 1 if self.deg_s == 1 else 6
        tab, tris, subarea, latx, dm = self._lattice_data(V, m)
        lv = vals[dm] @ tab.T  # (ncells, npts)
        sv = lv[:, tris]  # (ncells, ntri, 3)
        frac = tri_area_superlevel(sv.reshape(-1, 3), level).reshape(sv.shape[:2])
        area = float(np.sum(frac * subarea[:, None]))
        mask = lv >= level
        xr = float(np.max(latx[mask])) if np.any(mask) else 0.0
        return self.comm.allreduce(area, op=MPI.SUM), self.comm.allreduce(xr, op=MPI.MAX)

    def scalar(self, form):
        return self.comm.allreduce(fem.assemble_scalar(form), op=MPI.SUM)

    def q_area(self, mask):
        self.ind.x.array[:] = mask.astype(float)
        return self.scalar(self.form_ind)


def qpoint_xy(S):
    """Physical coordinates of the history (quadrature) points, ordered like S.Q dofs (affine cells)."""
    msh = S.msh
    X = S.Q.element.interpolation_points
    gx = msh.geometry.x[:, :2]
    gdm = msh.geometry.dofmap
    X0, X1, X2 = gx[gdm[:, 0]], gx[gdm[:, 1]], gx[gdm[:, 2]]
    P = X0[:, None, :] + X[None, :, 0:1] * (X1 - X0)[:, None, :] + X[None, :, 1:2] * (X2 - X0)[:, None, :]
    dm = S.Q.dofmap.list
    xy = np.zeros((len(S.H.x.array), 2))
    xy[dm.ravel()] = P.reshape(-1, 2)
    return xy


# ----------------------------------------------------------------------------------------
# Export helpers (serial)
# ----------------------------------------------------------------------------------------
def gather_cells(V, vals, static=False):
    """Gather per-cell data of a (collapsed) space on rank 0 (owned cells only).
    static=True: cell vertex coordinates and dof coordinates; else the per-cell dof values
    (shape (M, nd) or (M, nd, bs))."""
    msh = V.mesh
    comm = msh.comm
    nloc = msh.topology.index_map(2).size_local
    dm = V.dofmap.list[:nloc]
    bs = V.dofmap.index_map_bs
    if static:
        gdm = msh.geometry.dofmap[:nloc]
        loc = dict(cell_vertices=msh.geometry.x[gdm][:, :, :2].astype(np.float64))
        if "Quadrature" not in type(V.ufl_element()).__name__:
            loc["cell_dof_xy"] = V.tabulate_dof_coordinates()[:, :2][dm]
        else:  # quadrature space: physical point coordinates
            X = V.element.interpolation_points
            v = loc["cell_vertices"]
            loc["cell_dof_xy"] = (v[:, None, 0, :] + X[None, :, 0:1] * (v[:, None, 1, :] - v[:, None, 0, :])
                                  + X[None, :, 1:2] * (v[:, None, 2, :] - v[:, None, 0, :]))
        g = {k: comm.gather(a, root=0) for k, a in loc.items()}
        if comm.rank == 0:
            out = {k: np.concatenate(a) for k, a in g.items()}
            out["degree"] = np.array(V.ufl_element().embedded_superdegree)
            return out
        return None
    a = vals.reshape(-1, bs)[dm] if bs > 1 else vals[dm]
    g = comm.gather(a, root=0)
    return np.concatenate(g) if comm.rank == 0 else None
