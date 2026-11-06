"""The smooth practice problem D0 without a crack (element-test reference).

D0 : Omega = [0,1]^2, no notch (grad d . n = 0 on the whole boundary, no d Dirichlet).
     c = c_b(t) on the top edge y = 1 only, zero flux on the other edges.
     u_x = u_y = 0 at (0, 0); u_y = 0 at (1, 0); otherwise traction-free.
Initial: u = 0, c = c0, d = 0, H = alpha = W = 0.  Same c_b(t) = c0 + A (1 - cos(2 pi t/T))/2 as D1.

  python d0.py --problem D0 --n 64 --form F1 --div 40 --ncyc 20 --Gc .. --alphaT .. --A .. --T .. --chi .. --out DIR
  --n: N (structured N x N squares, crossed diagonals)
"""
import argparse
import os

import numpy as np

from solver import Params
import problem as pb
from problem import Problem

HOST_N = (8, 16, 32, 64)


def make_mesh(n, diag="crossed"):
    return pb.structured_rectangle((0.0, 0.0), (1.0, 1.0), int(n), int(n), diag)


COMMON = dict(make_mesh=make_mesh, notch_segments=[], tips=[], grid_bbox=((0.0, 1.0), (0.0, 1.0)), grid_n=101,
              near_dist=0.25, host_structured=HOST_N, host_series=True)

D0 = Problem(name="d0", u_points=[((0.0, 0.0), (0, 1)), ((1.0, 0.0), (1,))],
             c_boundary=lambda x: np.isclose(x[1], 1.0),
             probes=[(0.5, 1.0), (0.0, 1.0), (1.0, 1.0), (0.5, 0.9), (0.5, 0.75), (0.5, 0.5), (0.5, 0.0),
                     (0.25, 1.0), (0.0, 0.75)],
             **COMMON)

PROBLEMS = {"D0": D0}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--problem", default="D0")
    ap.add_argument("--n", default="64")
    ap.add_argument("--form", default="F1")
    ap.add_argument("--div", type=int, default=40)
    ap.add_argument("--ncyc", type=int, default=20)
    ap.add_argument("--Gc", type=float, required=True)
    ap.add_argument("--alphaT", type=float, required=True)
    ap.add_argument("--A", type=float, required=True)
    ap.add_argument("--T", type=float, required=True)
    ap.add_argument("--chi", type=float, default=0.25)
    ap.add_argument("--out", required=True)
    ap.add_argument("--export_dofs", action="store_true")
    ap.add_argument("--save_hist", action="store_true")
    ap.add_argument("--no_host_series", action="store_true")
    a = ap.parse_args()
    P = PROBLEMS[a.problem]
    if a.no_host_series:
        P.host_series = False
    prm = Params(Gc=a.Gc, alphaT=a.alphaT, A=a.A, T=a.T, chi=a.chi)
    pb.run(P, int(a.n), a.form, a.div, a.ncyc, prm, a.out, export_dofs=a.export_dofs,
           tag=os.path.basename(a.out.rstrip("/")), save_hist=a.save_hist)


if __name__ == "__main__":
    main()
