"""Record source ownership, generator transformations and environment versions."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import numpy
import scipy
import sympy
import abaqus_ufl

ROOT=Path(__file__).resolve().parents[1]


if __name__=="__main__":
    c=json.loads((ROOT/"element/uel_contract.json").read_text())
    files=c["source"]
    lines={s:len((ROOT/s).read_text().splitlines()) for s in files}
    # Seven inserted lines (including their separating blank lines) and one
    # replacement small-strain annotation are handwritten inside the kernel.
    injected=8
    generated=lines["element/generated_kernel.for"]-injected
    python_sources=[ROOT/"element/model.py",ROOT/"element/build.py"]+sorted((ROOT/"verification").glob("*.py"))
    assert Path(abaqus_ufl.__file__).resolve().is_relative_to(ROOT/"tools/abaqus_ufl")
    result=dict(
        recorded_at=subprocess.check_output(["date","-u","+%Y-%m-%dT%H:%M:%SZ"],text=True).strip(),
        python=sys.version,python_executable=sys.executable,
        abaqus_ufl_import_path=abaqus_ufl.__file__,
        numpy=numpy.__version__,scipy=scipy.__version__,sympy=sympy.__version__,
        gfortran=subprocess.check_output(["gfortran","--version"],text=True).splitlines()[0],
        generator="bundled abaqus_ufl, organizer-pinned public commit a02ed7c",
        generator_package_modified=False,host_modified=False,
        generated_lines=generated,handwritten_lines=sum(lines.values())-generated,
        handwritten_lines_inside_kernel=injected,
        generated_fraction=generated/sum(lines.values()),line_counts=lines,
        counting_convention="Physical source lines including comments and blank lines. Generated source includes unused shipped helper routines. The Python equations and formulation decisions were handwritten.",
        transformations=["Rename UEL entry and end to CHEMO_KERNEL",
            "Remove det(I+grad u)>0 finite-strain guard: fixed small-strain model has no such algebraic domain restriction",
            "Insert local-mode declarations, solve, strain correction and exact condensed Schur term"],
        handwritten_contributions=["Standard signature, validation, DTIME=0 and unsupported-request guards, transactional saved state",
            "Four zero-average element-local strain modes and pivoted static condensation, including every current-field coupling",
            "27-line row-sum lumped damage reaction driven by start-of-increment H; exact consistent diagonal tangent"],
        source_sha256={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in files},
        contract_sha256=hashlib.sha256((ROOT/"element/uel_contract.json").read_bytes()).hexdigest(),
        python_source_sha256={str(s.relative_to(ROOT)):hashlib.sha256(s.read_bytes()).hexdigest() for s in python_sources},
        abi="Standard Abaqus UEL signature; GNU Fortran/f2py through the unmodified host; no Abaqus execution",
        single_agent=True,network_used=False,package_installs=False)
    (ROOT/"verification/provenance.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
