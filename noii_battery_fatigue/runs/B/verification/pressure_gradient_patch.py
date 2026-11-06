"""Exact affine-rectangle pressure-gradient patch from u_x=a*x*y.

This independent direct compiled test exercises the displacement-gradient
contribution to pressure, unlike the clamped-compression drift oracle.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from common import ROOT, direct
from psilag_phase_oracle import compiled


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--element-root', type=Path, default=ROOT)
    parser.add_argument('--out', type=Path, default=ROOT/'verification/pressure_gradient_patch.json')
    args = parser.parse_args()
    root = args.element_root.resolve()
    assert root == ROOT or ROOT in root.parents
    module, contract = compiled(root)
    names = contract["dof_names"]
    dpn = len(names)
    index = {name: names.index(name) for name in names}
    p = dict(zip(contract["prop_names"], contract["props"]))
    Lx, Ly, amplitude, damage = 1.3, .8, .03, .27
    X = np.array([[0., 0.], [Lx, 0.], [Lx, Ly], [0., Ly]])
    mu = p["E"]/(2*(1+p["nu"]))
    lam = p["E"]*p["nu"]/((1+p["nu"])*(1-2*p["nu"]))
    pressure_coefficient = 2*(mu+lam)/3
    degradation = p["kappa"]+(1-p["kappa"])*(1-damage)**2
    q_slope = degradation*pressure_coefficient*amplitude
    U = np.zeros((4, dpn))
    U[:, index["ux"]] = amplitude*X[:, 0]*X[:, 1]
    U[:, index["c"]] = p["c0"]
    U[:, index["d"]] = damage
    U[:, index["q"]] = q_slope*X[:, 1]
    intact = next((name for name in ("vfrac", "v") if name in index), None)
    if intact is not None:
        U[:, index[intact]] = 1-damage
    rhs, tangent, states, cutback = direct(module, contract, X, U.ravel())
    if cutback != 1:
        raise AssertionError("valid pressure-gradient patch was rejected")
    R = -rhs.reshape(4, dpn)
    Jy = -p["chi"]*p["c0"]*q_slope
    expected_c = Jy*Lx/2*np.array([-1., -1., 1., 1.])
    pressure_residual = float(np.max(np.abs(R[:, index["q"]])))
    concentration_error = float(np.max(np.abs(R[:, index["c"]]-expected_c)))
    # Closed-form integral N_a*grad(N_b) on a Cartesian rectangle.
    sx = np.array([-1., 1., 1., -1.])
    sy = np.array([-1., -1., 1., 1.])
    derivative_error = 0.0
    for a in range(4):
        for b in range(4):
            Ix = sx[b]*Ly/8*(1+sy[a]*sy[b]/3)
            Iy = sy[b]*Lx/8*(1+sx[a]*sx[b]/3)
            row = a*dpn+index["q"]
            for component, integral in (("ux", Ix), ("uy", Iy)):
                expected = -degradation*pressure_coefficient*integral
                actual = tangent[row, b*dpn+index[component]]
                derivative_error = max(derivative_error, abs(actual-expected))
    assert pressure_residual < 1e-12, pressure_residual
    assert concentration_error < 1e-12, concentration_error
    assert derivative_error < 1e-11, derivative_error
    assert np.isfinite(tangent).all() and np.isfinite(states).all()
    # Remove the displacement gradient while retaining the same nonzero q.
    # An implementation that never computes div(u) would miss this change.
    broken = U.copy()
    broken[:, index["ux"]] = 0.
    broken_rhs = direct(module, contract, X, broken.ravel())[0].reshape(4, dpn)
    broken_pressure = float(np.max(np.abs(broken_rhs[:, index["q"]])))
    assert broken_pressure > 1e-3, broken_pressure
    result = dict(
        scope="exact displacement-gradient pressure and drift patch",
        dof_names=names, Lx=Lx, Ly=Ly, amplitude=amplitude, damage=damage,
        q_gradient=[0., q_slope], concentration_flux=[0., Jy],
        expected_c_residual=expected_c.tolist(),
        actual_c_residual=R[:, index["c"]].tolist(),
        q_residual_inf=pressure_residual,
        c_residual_error_inf=concentration_error,
        q_u_tangent_error_inf=derivative_error,
        removed_displacement_broken_control_q_residual=broken_pressure,
        source_sha256={name:hashlib.sha256((root/name).read_bytes()).hexdigest()
                       for name in contract["source"]+['element/uel_contract.json']},
        intact_constraint_inf=(None if intact is None else float(np.max(
            np.abs((U[:, index["d"]]-1)+U[:, index[intact]])))))
    args.out.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
