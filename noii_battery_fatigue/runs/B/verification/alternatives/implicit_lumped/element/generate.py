"""Generate the unmodified-tool UEL, and rename its entry for a checked wrapper.

Run from workspace root with PATH, PYTHONPATH, TMPDIR set as in README.
No changes are made to the abaqus_ufl package or its templates.
"""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parents[1]
sys.path[:0] = [str(WORKSPACE / "tools/abaqus_ufl"), str(ROOT / "element")]
from abaqus_ufl import generate_uel
from model import ChemoProblem, verification_states


def main():
    problem = ChemoProblem()
    checks = {}
    for label, state in verification_states():
        print(f"Verifying {label}")
        checks[label] = bool(problem.verify(state=state, verbose=True))
        if not checks[label]:
            raise RuntimeError(f"Verification failed: {label}")
    target = ROOT / "element/generated_model.for"
    generate_uel(problem, str(target), element="quad4", formulation="standard", mat_prefix="chem",
                 verify=dict(state=next(verification_states())[1]))
    src = target.read_text()
    assert src.count("SUBROUTINE UEL(") == 1
    src = src.replace("SUBROUTINE UEL(", "SUBROUTINE CHEM_UEL_CORE(", 1)
    # Rename the corresponding generated END SUBROUTINE as well.
    src = src.replace("END SUBROUTINE UEL", "END SUBROUTINE CHEM_UEL_CORE")
    from lump_reaction import transform
    raw_generated_lines = len(src.splitlines())
    src, extension = transform(src)
    target.write_text(src)
    manifest = dict(generator="abaqus_ufl", pinned_commit="a02ed7c",
                    package_modifications=[],
                    output_transformations=["Rename UEL entry to CHEM_UEL_CORE for checked wrapper", extension],
                    python_tangent_checks=checks,
                    sha256=hashlib.sha256(target.read_bytes()).hexdigest(),
                    generated_lines_before_output_extension=raw_generated_lines,
                    postprocessed_total_lines=len(src.splitlines()),
                    assembly_extension_original_lines_replaced=5,
                    assembly_extension_replacement_lines=10,
                    assembly_extension_added_lines=len(src.splitlines())-raw_generated_lines-5,
                    handwritten_wrapper_lines=len((ROOT / "element/uel.f90").read_text().splitlines()))
    (ROOT / "verification/generation.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
