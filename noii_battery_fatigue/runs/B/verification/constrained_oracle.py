"""Independent stationary compression/drift oracle on the unmodified host.

Run from the workspace with its pinned Python environment.  This script does
not import element/model.py or the generator.  u is constrained at every node,
c is prescribed only at left/right, and d and the projected q remain free.
A single DTIME=1e12 step approximates the steady problem; an exact-profile
residual check uses DU=0, hence has identically zero storage.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "host"))
import fehost  # noqa: E402
import meshes  # noqa: E402


def exact_profile(x, parameters, c_left, c_right):
    """F(c)=D*c+chi*B*c**2/2 is affine in x, B=4*beta*(mu+lam)/3."""
    E, nu = parameters["E"], parameters["nu"]
    mu, lam = E / (2 * (1 + nu)), E * nu / ((1 + nu) * (1 - 2 * nu))
    B = 4 * parameters["beta"] * (mu + lam) / 3
    a, D = parameters["chi"] * B, parameters["D"]
    F_left, F_right = D * c_left + a * c_left**2 / 2, D * c_right + a * c_right**2 / 2
    F = F_left + (F_right - F_left) * np.asarray(x)
    # The rationalized positive root avoids cancellation when chi is small.
    c = 2 * F / (D + np.sqrt(D**2 + 2 * a * F))
    dc_dx = (F_right - F_left) / (D + a * c)
    q = -B * (c - parameters["c0"])
    return c, q, dc_dx, F_right - F_left, B


def global_residual(element, nodes, elems, U):
    dpn = element.dof_per_node
    gm = (elems[:, :, None] * dpn + np.arange(dpn)).reshape(len(elems), -1)
    re, _ = element.element_rk_batch(nodes[elems], U[gm], np.zeros_like(U[gm]))
    assembled = np.zeros_like(U)
    np.add.at(assembled, gm.ravel(), re.ravel())
    return assembled


def independent_error_integrals(nodes, elems, U, element, parameters, c_left, c_right):
    """Geometry/shape/flux calculations written independently of the UEL."""
    uu = U.reshape(len(nodes), element.dof_per_node)
    ci, qi = element.comp("c"), element.comp("q")
    points, weights = np.polynomial.legendre.leggauss(4)
    totals = dict(c_L2=0.0, grad_c_L2=0.0, q_L2=0.0,
                  grad_q_L2=0.0, flux_L2=0.0, omitted_drift_flux_L2=0.0)
    for cell in elems:
        xy, cv, qv = nodes[cell], uu[cell, ci], uu[cell, qi]
        for i, xi in enumerate(points):
            for j, eta in enumerate(points):
                N = np.array([(1-xi)*(1-eta), (1+xi)*(1-eta),
                              (1+xi)*(1+eta), (1-xi)*(1+eta)]) / 4
                derivatives = np.array([[-(1-eta), -(1-xi)],
                                        [1-eta, -(1+xi)],
                                        [1+eta, 1+xi],
                                        [-(1+eta), 1-xi]]) / 4
                jacobian = xy.T @ derivatives
                det = np.linalg.det(jacobian)
                if det <= 0:
                    raise AssertionError("oracle encountered invalid geometry")
                gradient = derivatives @ np.linalg.inv(jacobian)
                x = float(N @ xy[:, 0])
                c, q, dc, flux, B = exact_profile(x, parameters, c_left, c_right)
                ch, qh = N @ cv, N @ qv
                gc, gq = cv @ gradient, qv @ gradient
                physical_flux = parameters["D"] * gc - parameters["chi"] * ch * gq
                exact_flux = np.array([flux, 0.0])
                differences = dict(
                    c_L2=(ch-c)**2,
                    grad_c_L2=np.sum((gc-np.array([dc, 0]))**2),
                    q_L2=(qh-q)**2,
                    grad_q_L2=np.sum((gq+np.array([B*dc, 0]))**2),
                    flux_L2=np.sum((physical_flux-exact_flux)**2),
                    omitted_drift_flux_L2=np.sum((parameters["D"]*gc-exact_flux)**2),
                )
                for name, value in differences.items():
                    totals[name] += weights[i] * weights[j] * det * value
    return {name: math.sqrt(value) for name, value in totals.items()}


def run_case(module, contract, mesh_name, n, seed, dt, c_left, c_right, tol,
             element_overrides=None):
    started = time.perf_counter()
    parameters = dict(zip(contract["prop_names"], contract["props"]))
    nodes, elems = (meshes.structured(n) if mesh_name == "structured"
                    else meshes.distorted(n, seed=seed))
    element = fehost.ContractUEL(module, contract, len(elems), dt)
    # Used only by the deliberately broken control: retain the original
    # analytical oracle while removing drift from the element under test.
    for name, value in (element_overrides or {}).items():
        element.props[contract["prop_names"].index(name)] = value
    dpn, ci, di, qi = element.dof_per_node, element.comp("c"), element.comp("d"), element.comp("q")
    c_exact, q_exact, _, flux_exact, B = exact_profile(nodes[:, 0], parameters, c_left, c_right)
    left, right = np.isclose(nodes[:, 0], 0), np.isclose(nodes[:, 0], 1)
    concentration_free = ~(left | right)
    exact_U = np.zeros((len(nodes), dpn))
    exact_U[:, ci], exact_U[:, qi] = c_exact, q_exact
    # The optional sixth field is the independently represented intact
    # fraction. Its constraint makes v=1 when this exact solution has d=0.
    intact_name = next((name for name in ("vfrac", "v")
                        if name in element.dof_names), None)
    if intact_name is not None:
        exact_U[:, element.comp(intact_name)] = 1.0
    exact_R = global_residual(element, nodes, elems, exact_U.ravel()).reshape(len(nodes), dpn)
    patch_c_residual = float(np.linalg.norm(exact_R[concentration_free, ci]))

    bc = {dpn * i + element.comp(component): 0.0
          for i in range(len(nodes)) for component in ("ux", "uy")}
    bc.update({dpn*i+ci: c_left for i in np.flatnonzero(left)})
    bc.update({dpn*i+ci: c_right for i in np.flatnonzero(right)})
    initial = np.zeros((len(nodes), dpn))
    initial[:, ci] = c_left + (c_right-c_left) * nodes[:, 0]
    initial[:, qi] = -B * (initial[:, ci]-parameters["c0"])
    solved, report = fehost.solve(
        nodes, elems, element, lambda fraction: bc, 1, dt,
        U0=initial.ravel(), tol=tol,
        residual_scale=fehost.residual_weights(element, contract))
    values = solved.reshape(len(nodes), dpn)
    residual = global_residual(element, nodes, elems, solved).reshape(len(nodes), dpn)
    errors = independent_error_integrals(nodes, elems, solved, element, parameters, c_left, c_right)
    node_error = float(np.max(np.abs(values[:, ci]-c_exact)))
    q_constraint = float(np.max(np.abs(values[:, qi]+B*(values[:, ci]-parameters["c0"]))))
    max_damage = float(np.max(np.abs(values[:, di])))
    v_constraint = None
    degradation_parity = None
    if intact_name is not None:
        constraint_values = (values[:, di]-1)+values[:, element.comp(intact_name)]
        v_constraint = float(np.max(np.abs(constraint_values)))
        direct_g = parameters["kappa"]+(1-parameters["kappa"])*(1-values[:, di])**2
        mixed_g = parameters["kappa"]+(1-parameters["kappa"])*values[:, element.comp(intact_name)]**2
        degradation_parity = float(np.max(np.abs(direct_g-mixed_g)))
    # An omitted-drift implementation must be distinguishable from this oracle.
    linear_profile_separation = float(np.max(np.abs(initial[:, ci]-c_exact)))
    if linear_profile_separation < 1e-3:
        raise AssertionError("chosen oracle is too weak to distinguish omitted drift")
    if not np.all(np.isfinite(solved)):
        raise AssertionError("non-finite solution")
    if q_constraint > 1e-7 or max_damage > 1e-8:
        raise AssertionError(f"oracle branch or q constraint failed: {q_constraint=}, {max_damage=}")
    if v_constraint is not None and v_constraint > 1e-8:
        raise AssertionError(f"intact-fraction constraint failed: {v_constraint=}")
    if mesh_name == "structured" and (node_error > 1e-6 or patch_c_residual > 1e-9):
        raise AssertionError(f"structured closed-form oracle failed: {node_error=}, {patch_c_residual=}")
    row = dict(mesh=mesh_name, n=n, seed=seed,
               exact_profile_free_c_residual=patch_c_residual,
               nodal_c_max_error=node_error, nodal_q_constraint_error=q_constraint,
               max_abs_d=max_damage, flux_exact=float(flux_exact),
               nodal_v_constraint_error=v_constraint,
               nodal_degradation_parity_error=degradation_parity,
               flux_left=float(-residual[left, ci].sum()),
               flux_right=float(residual[right, ci].sum()),
               omitted_drift_nodal_separation=linear_profile_separation,
               wall_seconds=time.perf_counter()-started, **errors,
               solver=report.summary())
    return row


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, nargs="+", default=[8, 16, 32])
    parser.add_argument("--element-root", type=Path, default=ROOT)
    parser.add_argument("--mesh", nargs="+", default=["structured", "distorted"],
                        choices=["structured", "distorted"])
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--dt", type=float, default=1e12)
    parser.add_argument("--tol", type=float, default=1e-11)
    parser.add_argument("--module-dir", type=Path)
    parser.add_argument("--module-name")
    parser.add_argument("--out", type=Path, default=ROOT / "verification" / "constrained_oracle_results")
    args = parser.parse_args(argv)
    root = args.element_root.resolve()
    assert root == ROOT or ROOT in root.parents
    contract = fehost.load_contract(str(root))
    if "q" not in contract["dof_names"]:
        raise ValueError("this oracle expects the proposed auxiliary q field")
    if (args.module_dir is None) != (args.module_name is None):
        raise ValueError("pass --module-dir and --module-name together")
    if args.module_dir is not None:
        sys.path.insert(0, str(args.module_dir.resolve()))
        module = importlib.import_module(args.module_name)
    else:
        build_dir = ROOT / "tmp" / "constrained_oracle_build" / str(os.getpid())
        module = fehost.build_module(str(root), contract, str(build_dir), f"uel_constrained_{os.getpid()}")
    parameters = dict(zip(contract["prop_names"], contract["props"]))
    c_left, c_right = parameters["c0"]+0.1, parameters["c0"]+0.8
    rows = []
    args.out.mkdir(parents=True, exist_ok=True)
    for mesh_name in args.mesh:
        for n in args.n:
            row = run_case(module, contract, mesh_name, n, args.seed,
                           args.dt, c_left, c_right, args.tol)
            rows.append(row)
            print(json.dumps(row), flush=True)
            (args.out / "results.json").write_text(json.dumps(dict(
                scope="independent constrained-compression steady drift oracle",
                dt=args.dt, c_left=c_left, c_right=c_right,
                parameters=parameters, results=rows), indent=2)+"\n")
    # Field error convergence should hold even though the structured nodal c
    # values are exact; those interpolants still have nonzero field errors.
    for mesh_name in args.mesh:
        subset = sorted((r for r in rows if r["mesh"] == mesh_name), key=lambda r: r["n"])
        for before, after in zip(subset, subset[1:]):
            if after["c_L2"] >= before["c_L2"] or after["flux_L2"] >= before["flux_L2"]:
                raise AssertionError(f"independent refinement failed for {mesh_name}")
    try:
        run_case(module, contract, "structured", 8, args.seed,
                 args.dt, c_left, c_right, args.tol,
                 element_overrides={"chi": 0.0})
    except AssertionError as failure:
        if "structured closed-form oracle failed" not in str(failure):
            raise
        broken_control = dict(omitted_drift="detected", detail=str(failure))
        print(json.dumps(dict(broken_control=broken_control)), flush=True)
    else:
        raise AssertionError("the deliberately omitted-drift control passed")
    (args.out / "results.json").write_text(json.dumps(dict(
        scope="independent constrained-compression steady drift oracle",
        dt=args.dt, c_left=c_left, c_right=c_right,
        parameters=parameters, results=rows,
        broken_control=broken_control,
        source_sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                       for name in contract["source"]+['element/uel_contract.json']}), indent=2)+"\n")
    scalar_names = [name for name in rows[0] if name != "solver"]
    with (args.out / "results.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=scalar_names)
        writer.writeheader()
        writer.writerows({key: row[key] for key in scalar_names} for row in rows)


if __name__ == "__main__":
    main()
