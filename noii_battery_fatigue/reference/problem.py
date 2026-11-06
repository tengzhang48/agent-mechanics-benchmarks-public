"""Problem definition, observer and time loop, independent of the geometry.

A Problem gives: a mesh builder, the pre-notch segments (d = 1 at every d-DOF on them), the
displacement point constraints, the facets with c = c_b(t) (default: the whole outer boundary),
sampling points (grid, probes, a ray from each notch tip for the crack front) and optional host-mesh
node sets.  d1.py defines D1; a new geometry needs only a new Problem (see README, "Reuse").

    from problem import Problem, run
    rows, snaps, S = run(problem, mesh_spec, form="F1", div=40, ncyc=20, prm=Params(...), out=DIR)
"""
import csv
import json
import os
import time
from dataclasses import dataclass, field
from typing import Callable, Optional, Sequence

import numpy as np
from mpi4py import MPI
import ufl
from dolfinx import fem
from dolfinx import mesh as dmesh

from solver import CMFSolver, Params, PointEvaluator, gather_cells, qpoint_xy, tri_area_superlevel

HERE = os.path.dirname(os.path.abspath(__file__))
FORMS = {"F1": dict(deg_u=2, deg_s=1, qdeg=4), "F2": dict(deg_u=3, deg_s=2, qdeg=8)}


# ----------------------------------------------------------------------------------------
# meshes
# ----------------------------------------------------------------------------------------
def mesh_from_npz(path):
    """Affine triangle mesh from an npz with x (N,2) and cells (M,3) (e.g. written by a gmsh script)."""
    import basix.ufl
    comm = MPI.COMM_WORLD
    m = np.load(path)
    el = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    cells = m["cells"].astype(np.int64) if comm.rank == 0 else np.zeros((0, 3), np.int64)
    x = m["x"] if comm.rank == 0 else np.zeros((0, 2))
    return dmesh.create_mesh(comm, cells, el, x)


def structured_rectangle(p0, p1, nx, ny, diag="crossed"):
    dg = {"crossed": dmesh.DiagonalType.crossed, "right": dmesh.DiagonalType.right,
          "left": dmesh.DiagonalType.left}[diag]
    return dmesh.create_rectangle(MPI.COMM_WORLD, [list(p0), list(p1)], [nx, ny],
                                  dmesh.CellType.triangle, diagonal=dg)


# ----------------------------------------------------------------------------------------
# geometry helpers
# ----------------------------------------------------------------------------------------
def dist_to_segments(xy, segments):
    """Distance of points xy (n,2) to the union of segments [((x0,y0),(x1,y1)), ...]."""
    d = np.full(len(xy), np.inf)
    for a, b in segments:
        a, b = np.asarray(a, float), np.asarray(b, float)
        ab = b - a
        s = np.clip(((xy - a) @ ab) / (ab @ ab), 0.0, 1.0)
        d = np.minimum(d, np.linalg.norm(xy - (a + s[:, None] * ab), axis=1))
    return d


@dataclass
class Problem:
    name: str
    make_mesh: Callable  # spec -> dolfinx mesh
    notch_segments: Sequence  # [((x0,y0),(x1,y1)), ...]; d = 1 at every d-DOF on them
    tips: Sequence  # [((xt,yt),(dx,dy),length)]: ray from each notch tip for the crack front
    u_points: Sequence  # [((x,y), (components...)), ...]: zero displacement components at points
    c_boundary: Optional[Callable] = None  # facet marker x -> bool; None = whole outer boundary
    grid_bbox: tuple = ((0.0, 1.0), (0.0, 1.0))
    grid_n: int = 101
    probes: Sequence = ()
    near_dist: float = 0.25  # "near the notch" = within this distance of the notch segments
    host_structured: Sequence = ()  # n for structured n x n quad host meshes of grid_bbox (squares only)
    host_series: bool = False  # store u, c, d at all host nodes at every accepted step (host_nodes.npz)
    notch_tol: float = 1e-10

    def on_notch(self, x):
        return dist_to_segments(np.asarray(x[:2]).T, self.notch_segments) < self.notch_tol

    def notch_dist(self, xy):
        return dist_to_segments(np.asarray(xy)[:, :2], self.notch_segments)


def _n_owned(dofs, V):
    return MPI.COMM_WORLD.allreduce(int(np.sum(np.asarray(dofs[1]) < V.dofmap.index_map.size_local)))


def setup_bcs(S, prm, P):
    msh, W = S.msh, S.W
    tdim = msh.topology.dim
    msh.topology.create_connectivity(tdim - 1, tdim)
    bcs = []
    Vu_c = [W.sub(0).sub(i).collapse()[0] for i in (0, 1)]
    zero = [fem.Function(V) for V in Vu_c]
    for (px, py), comps in P.u_points:
        mk = (lambda px, py: (lambda x: np.isclose(x[0], px) & np.isclose(x[1], py)))(px, py)
        for i in comps:
            dofs = fem.locate_dofs_geometrical((W.sub(0).sub(i), Vu_c[i]), mk)
            assert _n_owned(dofs, Vu_c[i]) == 1, f"point {px, py} is not a mesh vertex"
            bcs.append(fem.dirichletbc(zero[i], dofs, W.sub(0).sub(i)))
    if P.c_boundary is None:
        facets = dmesh.exterior_facet_indices(msh.topology)
    else:
        facets = dmesh.locate_entities_boundary(msh, tdim - 1, P.c_boundary)
    cb = fem.Function(S.Vc)
    cdofs = fem.locate_dofs_topological((W.sub(1), S.Vc), tdim - 1, facets)
    bcs.append(fem.dirichletbc(cb, cdofs, W.sub(1)))
    one = fem.Function(S.Vd)
    one.x.array[:] = 1.0
    ddofs = fem.locate_dofs_geometrical((W.sub(2), S.Vd), P.on_notch)
    bcs.append(fem.dirichletbc(one, ddofs, W.sub(2)))
    S.n_notch_dofs = _n_owned(ddofs, S.Vd)
    S.notch_dofs_d = np.asarray(ddofs[1])  # local indices in the collapsed d space

    def upd(t):
        cb.x.array[:] = prm.c_b(t)

    S.bcs = bcs
    S.bc_updaters = [upd]


def initial_state(S, prm):
    a = S.w_old.x.array
    a[:] = 0.0
    a[S.map_c] = prm.c0
    a[S.map_d[S.notch_dofs_d]] = 1.0
    S.w.x.array[:] = a


def host_areas(dn, level, h, m=8):
    """Reference counterparts of A_damage measured on a structured quad host mesh (spacing h) from
    the reference d at its nodes: (a) lumped nodal count sum_i [d_i >= level] |cell share of node i|,
    (b) area of the superlevel set of the bilinear (Q1) interpolant (m x m sub-cells per quad)."""
    w = np.full(dn.shape, h * h)
    w[0, :] *= 0.5; w[-1, :] *= 0.5; w[:, 0] *= 0.5; w[:, -1] *= 0.5
    a_nodal = float(np.sum(w * (dn >= level)))
    s = np.linspace(0.0, 1.0, m + 1)
    S_, T_ = np.meshgrid(s, s, indexing="xy")
    v00, v10, v01, v11 = dn[:-1, :-1], dn[:-1, 1:], dn[1:, :-1], dn[1:, 1:]  # [iy, ix]
    vals = (v00[..., None, None] * (1 - S_) * (1 - T_) + v10[..., None, None] * S_ * (1 - T_)
            + v01[..., None, None] * (1 - S_) * T_ + v11[..., None, None] * S_ * T_)
    a, b = vals[..., :-1, :-1], vals[..., :-1, 1:]
    c, d = vals[..., 1:, :-1], vals[..., 1:, 1:]
    t1 = np.stack([a, b, d], -1).reshape(-1, 3)
    t2 = np.stack([a, d, c], -1).reshape(-1, 3)
    sub = h * h / (m * m) / 2
    a_q1 = float((tri_area_superlevel(t1, level).sum() + tri_area_superlevel(t2, level).sum()) * sub)
    return a_nodal, a_q1


def _pname(p):
    return f"d_at_{p[0]:.4f}_{p[1]:.4f}"


class Observer:
    """Observables at every accepted step + field snapshots at the end of every cycle."""

    def __init__(self, S, prm, P):
        self.S, self.prm, self.P = S, prm, P
        self.comm = S.comm
        off = np.ones(len(S.map_d), dtype=bool)
        off[S.notch_dofs_d] = False
        self.off = off
        (x0, x1), (y0, y1) = P.grid_bbox
        self.grid_x = np.linspace(x0, x1, P.grid_n)
        self.grid_y = np.linspace(y0, y1, P.grid_n)
        X, Y = np.meshgrid(self.grid_x, self.grid_y, indexing="xy")  # [iy, ix]
        self.grid_eval = PointEvaluator(S.msh, np.c_[X.ravel(), Y.ravel()], allow_outside=True)
        self.grid_shape = X.shape
        self.rays = []
        for (tx, ty), (dx, dy), L in P.tips:
            s = np.linspace(0.0, L, 1401)
            pts = np.c_[tx + s * dx, ty + s * dy]
            self.rays.append((s, pts, PointEvaluator(S.msh, pts)))
        self.host = {}
        for nh in P.host_structured:
            gx = np.linspace(x0, x1, nh + 1)
            gy = np.linspace(y0, y1, nh + 1)
            HX, HY = np.meshgrid(gx, gy, indexing="xy")
            Pp = np.c_[HX.ravel(), HY.ravel()]
            self.host[nh] = (PointEvaluator(S.msh, Pp), P.on_notch(Pp.T), (x1 - x0) / nh)
        self.host_hist = {nh: [] for nh in self.host} if P.host_series else None
        self.host_xy = {nh: np.c_[np.meshgrid(np.linspace(x0, x1, nh + 1), np.linspace(y0, y1, nh + 1),
                                              indexing="xy")[0].ravel(),
                                  np.meshgrid(np.linspace(x0, x1, nh + 1), np.linspace(y0, y1, nh + 1),
                                              indexing="xy")[1].ravel()] for nh in self.host}
        self.host_t = []
        self.probe_pts = np.array(P.probes, dtype=float).reshape(-1, 2)
        self.probe_eval = PointEvaluator(S.msh, self.probe_pts) if len(self.probe_pts) else None
        self.d0 = None
        one = fem.form(fem.Constant(S.msh, 1.0) * ufl.dx(domain=S.msh), jit_options=S.jit())
        self.area = S.scalar(one)
        self.C_prev = prm.c0 * self.area
        self.cum_flux = 0.0
        self.flux_prev = 0.0
        self.alpha_first = None
        self.alpha_first_near = None
        self.near = P.notch_dist(qpoint_xy(S)) <= P.near_dist
        self.near_d = P.notch_dist(S.Vd.tabulate_dof_coordinates()) <= P.near_dist

    def gmax(self, a):
        return self.comm.allreduce(float(np.max(a)) if len(a) else -np.inf, op=MPI.MAX)

    def gmin(self, a):
        return self.comm.allreduce(float(np.min(a)) if len(a) else np.inf, op=MPI.MIN)

    def __call__(self, step, t, dt, its):
        S, prm, P = self.S, self.prm, self.P
        S.sync_collapsed()
        dv = S.d_values()
        if self.d0 is None:
            self.d0 = dv.copy()  # step 1: d driven by H = 0 only (baseline profile of the notch)
        A_dam, x_reach = S.superlevel(dv, S.Vd, 0.5)
        D_int = S.scalar(S.form_D)
        C_int = S.scalar(S.form_C)
        E_el = S.scalar(S.form_E)
        flux = S.scalar(S.form_flux)  # int_dOmega J_c . n ds (inflow rate, direct evaluation)
        q_cons = (C_int - self.C_prev) / dt  # discrete inflow rate (= reaction of the c rows)
        self.cum_flux += 0.5 * dt * (flux + self.flux_prev)  # trapezoid of the direct flux
        self.flux_prev = flux
        self.C_prev = C_int
        al = S.alpha.x.array
        a_over = S.q_area(al > prm.alphaT)
        a_over_near = S.q_area((al > prm.alphaT) & self.near)
        if self.alpha_first is None and a_over > 0:
            self.alpha_first = t
        if self.alpha_first_near is None and a_over_near > 0:
            self.alpha_first_near = t
        row = dict(step=step, t=t, cycle=t / prm.T, c_b=prm.c_b(t),
                   d_max=self.gmax(dv[self.off]),
                   A_damage=A_dam, D_int=D_int, C_int=C_int, E_el=E_el,
                   flux_in_direct=flux, flux_in_consistent=q_cons,
                   cum_flux_direct=self.cum_flux, C_gain=C_int - prm.c0 * self.area,
                   d_min=self.gmin(dv), d_max_all=self.gmax(dv),
                   d_excess_max=self.gmax(dv - self.d0), d_excess_min=self.gmin(dv - self.d0),
                   d_excess_max_near=self.gmax((dv - self.d0)[self.near_d]),
                   d_excess_max_far=self.gmax((dv - self.d0)[~self.near_d]),
                   d_max_far=self.gmax(dv[~self.near_d]))
        # crack front along each tip ray: first point where d < 0.5 (front of the region connected
        # to the notch), and last ray point with d >= 0.5
        for i, (s, pts, ev) in enumerate(self.rays):
            dl = ev(S.df).ravel()
            below = np.nonzero(dl < 0.5)[0]
            if len(below) == 0:
                sf = s[-1]
            elif below[0] == 0:
                sf = 0.0
            else:
                k = below[0] - 1
                sf = s[k] + (dl[k] - 0.5) / (dl[k] - dl[k + 1]) * (s[k + 1] - s[k])
            above = np.nonzero(dl >= 0.5)[0]
            sfx = "" if i == 0 else f"_{i}"
            (tx, ty), (dx, dy), L = P.tips[i]
            row[f"front_s{sfx}"] = float(sf)
            row[f"x_tip{sfx}"] = float(tx + sf * dx)
            row[f"y_tip{sfx}"] = float(ty + sf * dy)
            row[f"last_s{sfx}"] = float(s[above[-1]]) if len(above) else 0.0
        row.update(x_reach=x_reach,
                   alpha_max=self.gmax(al), alpha_max_near=self.gmax(al[self.near]),
                   alpha_max_far=self.gmax(al[~self.near]),
                   area_alpha_gt_T=a_over, area_alpha_gt_T_near=a_over_near,
                   H_max=self.gmax(S.H.x.array), H_dmin=S.hist_dmin[0], alpha_dmin=S.hist_dmin[1],
                   c_min=self.gmin(S.cf.x.array), c_max=self.gmax(S.cf.x.array),
                   newton_its=its)
        if self.host_hist is not None:
            self.host_t.append(t)
        for nh, (ev, notch, h) in self.host.items():
            dh = ev(S.df).ravel()
            if self.host_hist is not None:
                uh = ev(S.uf)
                self.host_hist[nh].append(np.stack([uh[:, 0], uh[:, 1], ev(S.cf).ravel(), dh.copy()]).astype(np.float32))
            dh[notch] = 1.0
            row[f"d_max_host{nh}"] = float(np.max(dh[~notch]))
            an, aq = host_areas(dh.reshape(nh + 1, nh + 1), 0.5, h)
            row[f"A_damage_host{nh}_nodal"] = an
            row[f"A_damage_host{nh}_q1"] = aq
        if self.probe_eval is not None:
            for p, val in zip(self.probe_pts, self.probe_eval(S.df).ravel()):
                row[_pname(p)] = float(val)
        return row

    def snapshot(self):
        S = self.S
        S.sync_collapsed()
        u = self.grid_eval(S.uf)
        sh = self.grid_shape
        return dict(ux=u[:, 0].reshape(sh), uy=u[:, 1].reshape(sh),
                    c=self.grid_eval(S.cf).reshape(sh), d=self.grid_eval(S.df).reshape(sh),
                    q=self.grid_eval(S.qf).reshape(sh))


def run(P, mesh_spec, form="F1", div=40, ncyc=15, prm=None, out=None, export_dofs=False,
        verbose=True, stop_if=None, tag=None, save_hist=False):
    prm = prm or Params()
    comm = MPI.COMM_WORLD
    t_start = time.time()
    msh = P.make_mesh(mesh_spec)
    S = CMFSolver(msh, prm, prefix=f"{P.name}_{form}_{mesh_spec}_", **FORMS[form])
    setup_bcs(S, prm, P)
    S.finalize()
    initial_state(S, prm)
    obs = Observer(S, prm, P)
    if out and comm.rank == 0:
        os.makedirs(out, exist_ok=True)
    rows, snaps, dofsnaps = [], [], []
    dt = prm.T / div
    nsteps = div * ncyc
    t_setup = time.time() - t_start
    static = None
    info = dict(problem=P.name, mesh=str(mesh_spec), form=form, div=div, ncyc=ncyc)
    for k in range(1, nsteps + 1):
        t_old, t_new = (k - 1) * dt, k * dt
        sub = S.advance(t_old, t_new)
        its = sum(s[1] for s in sub)
        row = obs(k, t_new, dt, its)
        row["substeps"] = len(sub)
        row["wall"] = time.time() - t_start
        rows.append(row)
        if k % div == 0:
            snaps.append(obs.snapshot())
            if export_dofs:
                if static is None:
                    static = {nm: gather_cells(V, None, static=True)
                              for nm, V in (("u", S.Vu), ("c", S.Vc), ("d", S.Vd), ("q", S.Vq))}
                S.sync_collapsed()
                dofsnaps.append({nm: gather_cells(V, f.x.array) for nm, V, f in
                                 (("u", S.Vu, S.uf), ("c", S.Vc, S.cf), ("d", S.Vd, S.df), ("q", S.Vq, S.qf))})
            if verbose and comm.rank == 0:
                print(f"[{tag or ''} {mesh_spec} {form} div={div}] cycle {k // div}: d_max={row['d_max']:.4f} "
                      f"dexN={row['d_excess_max_near']:.4f} dexF={row['d_excess_max_far']:.4f} "
                      f"A={row['A_damage']:.5f} D={row['D_int']:.5f} C={row['C_int']:.5f} E={row['E_el']:.3e} "
                      f"xtip={row.get('x_tip', np.nan):.4f} xreach={row['x_reach']:.4f} "
                      f"amax_near={row['alpha_max_near']:.3e} amax_far={row['alpha_max_far']:.3e} "
                      f"Aa={row['area_alpha_gt_T']:.4f} Hmax={row['H_max']:.3e} its={its} "
                      f"wall={row['wall']:.0f}s", flush=True)
            if out:
                write_outputs(out, rows, snaps, dofsnaps, static, S, obs, prm, info, t_setup, save_hist)
        if stop_if is not None and stop_if(row):
            break
    if out:
        write_outputs(out, rows, snaps, dofsnaps, static, S, obs, prm, info, t_setup, save_hist)
    return rows, snaps, S


def write_outputs(out, rows, snaps, dofsnaps, static, S, obs, prm, info, t_setup, save_hist):
    comm = MPI.COMM_WORLD
    if save_hist:
        hist = {nm: gather_cells(S.Q, f.x.array) for nm, f in (("H", S.H), ("alpha", S.alpha), ("W", S.Wh))}
        hist_xy = gather_cells(S.Q, None, static=True)
    if comm.rank != 0:
        return
    keys = list(rows[0].keys())
    with open(os.path.join(out, "observables.csv"), "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=keys)
        wr.writeheader()
        for r in rows:
            wr.writerow(r)
    if snaps:
        np.savez_compressed(os.path.join(out, "fields_grid101.npz"),
                            x=obs.grid_x, y=obs.grid_y,
                            t=np.array([(i + 1) * prm.T for i in range(len(snaps))]),
                            **{k: np.array([s[k] for s in snaps]) for k in snaps[0]},
                            note="arrays [cycle, iy, ix] at point (x[ix], y[iy]); end of each cycle; "
                                 "NaN outside the domain")
    if dofsnaps:
        sp = {}
        for nm in ("u", "c", "d", "q"):
            for k, v in static[nm].items():
                sp[f"{nm}_{k}"] = v
            sp[f"{nm}_values"] = np.array([s[nm] for s in dofsnaps])
        sp["t"] = np.array([(i + 1) * prm.T for i in range(len(dofsnaps))])
        sp["note"] = ("per-cell Lagrange data: X_cell_vertices (M,3,2), X_cell_dof_xy (M,nd,2), "
                      "X_values (ncycle,M,nd[,2]); evaluate with evaluate_reference.py")
        np.savez_compressed(os.path.join(out, "fields_dofs.npz"), **sp)
    if obs.host_hist is not None and obs.host_t:
        hh = {f"n{nh}_values": np.array(v) for nh, v in obs.host_hist.items()}
        hh.update({f"n{nh}_xy": obs.host_xy[nh] for nh in obs.host_hist})
        np.savez_compressed(os.path.join(out, "host_nodes.npz"), t=np.array(obs.host_t), **hh,
                            note="nX_values[step, field, node] float32, fields (u_x, u_y, c, d); node k at "
                                 "nX_xy[k] = (i/n, j/n), k = j*(n+1) + i; every accepted step")
    if save_hist:
        np.savez_compressed(os.path.join(out, "history_end.npz"), xy=hist_xy["cell_dof_xy"],
                            **{k: v for k, v in hist.items()})
    log = [dict(t=e["t"], dt=e["dt"], its=e["its"], reason=e["reason"], res=e["res"],
                depth=e["depth"], failed=e.get("failed", False)) for e in S.newton_log]
    with open(os.path.join(out, "newton_log.json"), "w") as fh:
        json.dump(log, fh)
    meta = dict(info, params=prm.asdict(),
                deg_u=S.deg_u, deg_s=S.deg_s, qdeg=S.qdeg, qdeg_obs=S.qdeg_obs, mpi_size=comm.size,
                ndofs=int(S.W.dofmap.index_map.size_global * S.W.dofmap.index_map_bs),
                ncells=int(S.msh.topology.index_map(2).size_global),
                n_notch_dofs_d=int(S.n_notch_dofs), domain_area=obs.area,
                t_setup=t_setup, wall=rows[-1]["wall"],
                steps_done=len(rows), alpha_first_exceed_t=obs.alpha_first,
                alpha_first_exceed_near_t=obs.alpha_first_near,
                failed_attempts=sum(1 for e in S.newton_log if e.get("failed")),
                newton_its_total=int(sum(r["newton_its"] for r in rows)))
    with open(os.path.join(out, "meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
