"""Deterministically generate the UEL kernel and contract from model.py.

Run from workspace root with: source env.sh; python element/build.py
Transformations of generator output: rename UEL to CHEMO_KERNEL and remove
the finite-strain det(F)>0 guard, which is not a restriction of this fixed
small-strain model. The handwritten accepted_damage.f90 adds the row-sum
lumped reaction driven by accepted start-of-increment history.
"""
from pathlib import Path
import sys
import json
import hashlib
import argparse

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from element.model import ChemoFatigueQuad, LumpedCoreQuad, verification_states
import abaqus_ufl as au


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "element"))
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    # Check the full pointwise physics and the selected generated core.
    for declaration in (ChemoFatigueQuad,LumpedCoreQuad):
        problem=declaration()
        states=verification_states(problem)
        for i,state in enumerate(states):
            print(f"Verifying {declaration.__name__}, regime {i}")
            problem.verify(state=state,verbose=True,tol=2e-5)
    generated = out / "generated_kernel.for"
    au.generate_uel(problem, str(generated), element="quad4",
                    formulation="standard", mat_prefix="cf",
                    verify=dict(state=states[0], tol=2e-5))
    text = generated.read_text()
    assert text.count("SUBROUTINE UEL(") == 1
    text = text.replace("SUBROUTINE UEL(", "SUBROUTINE CHEMO_KERNEL(")
    # The generator's end statement is unnamed in the current pinned version.
    text = text.replace("END SUBROUTINE UEL", "END SUBROUTINE CHEMO_KERNEL")
    start=text.index("C       Reject invalid deformation states: det F must stay")
    end=text.index("        END IF",start)+len("        END IF")
    text=text[:start]+"C       Small-strain model: no det(I+grad u) restriction."+text[end:]
    generated.write_text(text)
    contract = dict(
        source=["element/generated_kernel.for", "element/accepted_damage.f90", "element/uel.f90"],
        compile_flags=["-ffixed-line-length-none"],
        dof_per_node=6, dof_names=["ux", "uy", "c", "d", "p", "s"],
        abaqus_dofs=[1, 2, 11, 12, 13, 14],
        prop_names=problem._mat.props_names, props=problem._mat.props_array.tolist(),
        jprops=[], nsvars=12, svars_initial=[0.0]*12,
        history_slots=dict(H=[1, 4, 7, 10], alpha=[2, 5, 8, 11]),
        lflags=[1, 0, 1, 0, 0, 0],
        # Keep mechanics and the pressure constraint resolved as g tends to
        # kappa; the driver uses weights equal to 1/these characteristic scales.
        residual_scale=dict(ux=0.01, uy=0.01, c=1.0, d=1.0, p=0.01, s=1.0),
        time_scheme="Backward Euler concentration, first-order accepted-history split",
        damage_reaction="Gauss-point H_n, row-sum lumped reaction on nodal intactness")
    (out / "uel_contract.json").write_text(json.dumps(contract, indent=2)+"\n")
    print(f"Generated {generated.relative_to(ROOT)}, {len(text.splitlines())} lines, "
          f"sha256={hashlib.sha256(text.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
