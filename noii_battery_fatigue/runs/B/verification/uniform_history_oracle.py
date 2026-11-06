"""Compiled UEL versus independent frozen uniform-strain history values.

The numerical expectations are hand/model-equation calculations documented in
independent_review.md. This file imports no generated model or assembly code.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from common import ROOT, direct, identical_bytes
from psilag_phase_oracle import compiled


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--element-root', type=Path, default=ROOT)
    parser.add_argument('--out', type=Path, default=ROOT/'verification/uniform_history_oracle.json')
    args = parser.parse_args()
    root = args.element_root.resolve()
    assert root == ROOT or ROOT in root.parents
    module, contract = compiled(root)
    names = contract["dof_names"]
    index = {name: names.index(name) for name in names}
    dpn = len(names)
    stride = contract['nsvars']//4
    assert stride in (3, 4)
    p = dict(zip(contract["prop_names"], contract["props"]))
    props = np.asarray(contract["props"], float).copy()
    props[contract["prop_names"].index("G_c")] = .002
    props[contract["prop_names"].index("alpha_T")] = .0005
    # Expected values are independent of element geometry: affine physical
    # displacement has constant strain even on this non-affine Q1 mapping.
    X = np.array([[0., 0.], [1.1, -.05], [.95, .85], [-.12, 1.03]])
    ee = np.array([[.03, .012], [.012, -.01]])
    concentration, damage = .55, .4
    g = .3600064
    exact_p = .11837121212121211

    def field(elastic_scale):
        U = np.zeros((4, dpn))
        strain = elastic_scale*ee+p["beta"]*(concentration-p["c0"])*np.eye(2)
        U[:, [index["ux"], index["uy"]]] = X@strain.T
        U[:, index["c"]] = concentration
        U[:, index["d"]] = damage
        U[:, index["q"]] = g*exact_p*elastic_scale
        intact = next((name for name in ("vfrac", "v") if name in index), None)
        if intact is not None:
            U[:, index[intact]] = 1-damage
        return U.ravel()

    expected = np.array([1.4492080900029878,
                         .0020754914424242423,
                         .0020754914424242423])
    expected_unload = expected.copy()
    expected_unload[2] *= .25
    expected_reload = np.array([2.1337751328249275,
                                .003632110024242424,
                                .0020754914424242423])
    expected_floor = np.array([6.918181818181818e14,
                               .10107549144242424,
                               .0020754914424242423])
    psi_hand = .005765151515151515

    def old_state(triple, psi):
        row = np.zeros(stride)
        row[:3] = triple
        if stride == 4:
            row[3] = psi
        return np.tile(row, 4)

    cases = [
        ("loading_hand_values", field(1), old_state([.3, .001, .001], psi_hand/2), expected, psi_hand),
        ("unloading_frozen_history", field(.5), old_state(expected, psi_hand), expected_unload, psi_hand/4),
        ("reloading_accumulated_increase", field(1), old_state(expected_unload, psi_hand/4), expected_reload, psi_hand),
        ("fatigue_floor", field(1), old_state([.3, .1, .001], psi_hand/2), expected_floor, psi_hand),
    ]
    results = {}
    for label, U, old, reference, expected_psi in cases:
        output = direct(module, contract, X, U, state=old, props=props)
        rhs, tangent, state, cutback = output
        assert cutback == 1, (label, cutback)
        expected_points = np.tile(reference, (4, 1))
        computed = state.reshape(4, stride)
        relative_error = float(np.max(abs(computed[:,:3]-expected_points)/np.maximum(1., abs(expected_points))))
        assert relative_error < 1e-12, (label, relative_error, computed, expected_points)
        assert np.isfinite(rhs).all() and np.isfinite(tangent).all()
        assert np.all(computed[:, :2] >= old.reshape(4, stride)[:, :2])
        if stride == 4:
            assert np.max(abs(computed[:,3]-expected_psi)) < 1e-13
        results[label] = dict(expected_H_alpha_W=reference.tolist(),
                             computed_first_point=computed[0].tolist(),
                             expected_current_psi=expected_psi if stride==4 else None,
                             max_point_scaled_error=relative_error)
    # A discarded high trial must not alter the exact same subsequent call.
    U, old = cases[0][1:3]
    baseline = direct(module, contract, X, U, state=old, props=props)
    for scale in (3., .1, 5., .25, 2.):
        direct(module, contract, X, field(scale), state=old, props=props)
    repeat = direct(module, contract, X, U, state=old, props=props)
    assert all(identical_bytes(a, b) for a, b in zip(baseline, repeat))
    results["discarded_trial_repeat"] = "bitwise-identical RHS, tangent, SVARS and PNEWDT"
    # A deliberately wrong history update accumulating full W on reload
    # would predict alpha=alpha_n+W, instead of alpha_n+(W-W_n).
    wrong_reload_alpha = expected[1]+expected[2]
    assert abs(wrong_reload_alpha-expected_reload[1]) > 5e-4
    erased_W_state = old_state(expected_unload, psi_hand/4).reshape(4,stride)
    erased_W_state[:, 2] = 0.
    erased_output = direct(module, contract, X, field(1),
                           state=erased_W_state.ravel(), props=props)[2].reshape(4, stride)
    assert np.max(np.abs(erased_output[:, 1]-expected_reload[1])) > 5e-4
    assert np.max(np.abs(erased_output[:, 1]-wrong_reload_alpha)) < 1e-14
    results["full_W_accumulation_broken_control"] = dict(
        wrong_alpha=wrong_reload_alpha, expected_alpha=expected_reload[1],
        separation=wrong_reload_alpha-expected_reload[1],
        compiled_erased_W_alpha=erased_output[0, 1])
    results["source_sha256"] = {
        name: hashlib.sha256((root/name).read_bytes()).hexdigest()
        for name in contract["source"]+['element/uel_contract.json']}
    args.out.write_text(json.dumps(results, indent=2)+"\n")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
