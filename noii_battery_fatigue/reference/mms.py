"""Manufactured-solution verification of the reference solver (F1 and F2).

Exact smooth fields u, c, d and a prescribed smooth history field H (no history update):
sources are built from the strong form with UFL differentiation of the exact expressions
(an independent code path from the weak residual), q_ex = g(d_ex) p_eff(u_ex, c_ex).
All of u, c, d are Dirichlet on the whole boundary; q has no boundary condition.

Two variants:
  steady    : time derivative dropped (time-independent problem)
  transient : c_ex = c_s(x) + t c_1(x) (linear in t) so backward Euler has no time
              truncation error; 4 steps of dt = 0.25 from the exact initial state.
Usage: python mms.py [F1|F2|all]
"""
import json
import os
import sys
import time

import numpy as np
import ufl
from mpi4py import MPI
from dolfinx import fem
from dolfinx import mesh as dmesh

from solver import CMFSolver, Params, kinematics

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", "mms")


def exact(x, t, prm):
    pi = np.pi
    ux = 0.03 * ufl.sin(1.2 * pi * x[0] + 0.3) * ufl.cos(0.8 * pi * x[1] - 0.2)
    uy = 0.02 * ufl.cos(0.9 * pi * x[0] + 0.1) * ufl.sin(1.1 * pi * x[1] + 0.4)
    u = ufl.as_vector([ux, uy])
    c = 0.45 + 0.15 * ufl.sin(2.1 * x[0] + 0.7) * ufl.cos(1.6 * x[1] - 0.3) \
        + t * 0.1 * ufl.cos(1.3 * x[0] - 0.5 * x[1])
    d = 0.35 + 0.25 * ufl.cos(1.9 * x[0] - 0.4) * ufl.sin(2.3 * x[1] + 0.5)
    H = 0.6 + 0.4 * ufl.sin(x[0] + 2.0 * x[1])
    c_t = 0.1 * ufl.cos(1.3 * x[0] - 0.5 * x[1])
    eps_e, sig_eff, g, p_eff = kinematics(u, c, d, prm)
    q = g * p_eff
    return u, c, d, q, H, c_t


def run_case(n, deg_u, deg_s, qdeg, transient, prm, diag=dmesh.DiagonalType.crossed):
    msh = dmesh.create_rectangle(MPI.COMM_WORLD, [[0.0, 0.0], [1.0, 1.0]], [n, n],
                                 dmesh.CellType.triangle, diagonal=diag)
    S = CMFSolver(msh, prm, deg_u=deg_u, deg_s=deg_s, qdeg=qdeg, steady=not transient,
                  prefix=f"mms{n}_{deg_u}{deg_s}{int(transient)}_", snes_rtol=1e-12, snes_atol=1e-14)
    x = ufl.SpatialCoordinate(msh)
    tc = fem.Constant(msh, 0.0)
    u_e, c_e, d_e, q_e, H_e, ct_e = exact(x, tc, prm)
    eps_e, sig_eff, g, p_eff = kinematics(u_e, c_e, d_e, prm)
    J_e = g * prm.D * ufl.grad(c_e) - prm.chi * c_e * ufl.grad(g * p_eff)
    f_u = -ufl.div(g * sig_eff)
    s_c = (ct_e if transient else 0.0) - ufl.div(J_e)
    s_d = -prm.ell ** 2 * ufl.div(ufl.grad(d_e)) + d_e - (1 - prm.kappa) * (1 - d_e) * H_e
    v, wc, wd, r = S.tests
    dxs = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": 12})
    scale = S.dt if transient else 1.0
    extra = -(ufl.inner(f_u, v) + scale * s_c * wc + s_d * wd) * dxs
    # prescribed history field at the quadrature points
    S.H.interpolate(fem.Expression(H_e, S.Q.element.interpolation_points))
    # Dirichlet data
    msh.topology.create_connectivity(1, 2)
    bf = dmesh.exterior_facet_indices(msh.topology)
    bcs, upd = [], []
    for sub, V, ex in [(0, S.Vu, u_e), (1, S.Vc, c_e), (2, S.Vd, d_e)]:
        gfun = fem.Function(V)
        expr = fem.Expression(ex, V.element.interpolation_points)
        dofs = fem.locate_dofs_topological((S.W.sub(sub), V), 1, bf)
        bcs.append(fem.dirichletbc(gfun, dofs, S.W.sub(sub)))
        upd.append((gfun, expr))

    def update(t):
        tc.value = t
        for gfun, expr in upd:
            gfun.interpolate(expr)

    S.bcs = bcs
    S.bc_updaters = [update]
    S.finalize(extra_F=extra)
    # initial state = interpolant of the exact fields at t = 0 (also initial guess)
    tc.value = 0.0
    for sub, V, ex in [(0, S.Vu, u_e), (1, S.Vc, c_e), (2, S.Vd, d_e), (3, S.Vq, q_e)]:
        f = fem.Function(V)
        f.interpolate(fem.Expression(ex, V.element.interpolation_points))
        idx = [S.map_u, S.map_c, S.map_d, S.map_q][sub]
        S.w_old.x.array[idx] = f.x.array
    t0 = time.time()
    its = []
    if transient:
        dt = 0.25
        for k in range(1, 5):
            ok, reason, it, res = S._solve(k * dt, dt)
            assert ok, (reason, res)
            S.w_old.x.array[:] = S.w.x.array  # no history update in MMS (H prescribed)
            its.append(it)
        tc.value = 1.0
    else:
        ok, reason, it, res = S._solve(0.0, 1.0)
        assert ok, (reason, res)
        its.append(it)
    wall = time.time() - t0
    uh, ch, dh, qh = ufl.split(S.w)
    dxe = ufl.Measure("dx", domain=msh, metadata={"quadrature_degree": 14})
    err = {}
    for name, fh, fe in [("u", uh, u_e), ("c", ch, c_e), ("d", dh, d_e), ("q", qh, q_e)]:
        e = fh - fe
        L2 = np.sqrt(fem.assemble_scalar(fem.form(ufl.inner(e, e) * dxe)))
        H1 = np.sqrt(fem.assemble_scalar(fem.form(ufl.inner(ufl.grad(e), ufl.grad(e)) * dxe)))
        nrm = np.sqrt(fem.assemble_scalar(fem.form(ufl.inner(fe, fe) * dxe)))
        err[name] = dict(L2=float(L2), H1=float(H1), norm=float(nrm))
    ndof = S.W.dofmap.index_map.size_global * S.W.dofmap.index_map_bs
    return dict(n=n, h=1.0 / n, ndof=int(ndof), its=its, wall=wall, err=err, res_last=res)


def rates(rows, field, norm):
    out = []
    for a, b in zip(rows[:-1], rows[1:]):
        out.append(np.log(a["err"][field][norm] / b["err"][field][norm]) / np.log(a["h"] / b["h"]))
    return out


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    os.makedirs(OUT, exist_ok=True)
    prm = Params()
    cases = []
    if which in ("F1", "all"):
        cases += [("F1", 2, 1, 4, False, [8, 16, 32, 64, 128]),
                  ("F1", 2, 1, 4, True, [8, 16, 32, 64, 128])]
    if which in ("P2P2", "all"):
        # diagnostic: u P2 with P2 scalar fields (shows the u L2 rate is limited by P1 c, d, q in F1)
        cases += [("P2P2", 2, 2, 6, False, [8, 16, 32, 64])]
    if which in ("F2", "all"):
        cases += [("F2", 3, 2, 8, False, [4, 8, 16, 32, 64]),
                  ("F2", 3, 2, 8, True, [4, 8, 16, 32, 64])]
    allres = {}
    for name, du, ds, qd, tr, ns in cases:
        key = f"{name}_{'transient' if tr else 'steady'}"
        rows = []
        for n in ns:
            r = run_case(n, du, ds, qd, tr, prm)
            rows.append(r)
            print(key, n, {k: (f"{v['L2']:.3e}", f"{v['H1']:.3e}") for k, v in r["err"].items()},
                  "its", r["its"], f"{r['wall']:.1f}s", flush=True)
        rt = {f: {nm: rates(rows, f, nm) for nm in ("L2", "H1")} for f in ("u", "c", "d", "q")}
        print(key, "rates", json.dumps(rt), flush=True)
        allres[key] = dict(rows=rows, rates=rt, deg_u=du, deg_s=ds, qdeg=qd)
        with open(os.path.join(OUT, f"mms_{which}.json"), "w") as fh:
            json.dump(allres, fh, indent=1)


if __name__ == "__main__":
    main()
