"""Record reproducible source/compiler identities without reading outside inputs."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def counts(path):
    lines=path.read_text().splitlines()
    fixed=path.suffix==".for"
    code=[x for x in lines if x.strip() and not (x[0] in "cC*!" if fixed else x.lstrip().startswith("!"))]
    return dict(total_lines=len(lines),nonblank_noncomment_lines=len(code))


def main():
    contract=json.loads((ROOT/"element/uel_contract.json").read_text())
    element={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/"element").glob("*")) if p.is_file()}
    host={str(p.relative_to(ROOT)):digest(p) for p in sorted((ROOT/"host").rglob("*"))
          if p.is_file() and p.suffix in (".py",".f90",".md",".json")}
    gen_root=ROOT/"tools/abaqus_ufl/abaqus_ufl"
    generator={str(p.relative_to(ROOT)):digest(p) for p in sorted(gen_root.rglob("*"))
               if p.is_file() and p.suffix in (".py",".for")}
    generated=counts(ROOT/"element/generated_model.for")
    handwritten=counts(ROOT/"element/uel.f90")
    generation=json.loads((ROOT/'verification/generation.json').read_text())
    assert generation['sha256']==element['element/generated_model.for']
    raw_lines=generation.get('generated_lines_before_output_extension',generated['total_lines'])
    replaced=generation.get('assembly_extension_original_lines_replaced',0)
    authored_replacements=generation.get('assembly_extension_replacement_lines',0)
    authored_additions=generation.get('assembly_extension_added_lines',0)
    # Conservatively classify both renamed entry/END lines as authored output
    # transformations, as well as every rewritten/added assembly line.
    renamed=2
    unchanged_generated=raw_lines-replaced-renamed
    authored_core=authored_replacements+authored_additions+renamed
    assert unchanged_generated+authored_core==generated['total_lines']
    delivered_total=generated['total_lines']+handwritten['total_lines']
    allocation=dict(unmodified_generated_lines=unchanged_generated,
                    rewritten_or_added_core_lines=authored_core,
                    handwritten_wrapper_lines=handwritten['total_lines'],
                    delivered_total_lines=delivered_total,
                    fraction_unmodified_generated=unchanged_generated/delivered_total,
                    counting_policy='Every modified core line, including entry renaming and comments, is classified as authored.')
    earlier_path=ROOT/'verification/implicit_active_source/verification/provenance.json'
    audit={}
    if earlier_path.exists():
        earlier=json.loads(earlier_path.read_text())
        audit=dict(comparator=str(earlier_path.relative_to(ROOT)),
                   host_identical=host==earlier['host_sha256'],
                   generator_identical=generator==earlier['generator_sha256'])
        assert audit['host_identical'] and audit['generator_identical']
    data=dict(python_executable=sys.executable,python_version=sys.version,
              gfortran_executable=shutil.which("gfortran"),
              gfortran_version=subprocess.run(["gfortran","--version"],capture_output=True,text=True,check=True).stdout.splitlines()[0],
              generator_pinned_revision_as_declared_by_organizer="a02ed7c",
              revision_independently_verified=False,revision_note="vendored files have no Git metadata",
              package_modifications=[],host_modifications=[],
              generator_output_transformations=generation['output_transformations'],
              postprocessed_fortran_core=generated,handwritten_fortran_wrapper=handwritten,
              fortran_line_allocation=allocation,
              generated_fraction_of_delivered_fortran_lines=allocation['fraction_unmodified_generated'],
              line_count_caveat="Generic unused shipped templates are included; percentage does not measure authorship of physics.",
              generated_responsibilities=["constitutive law translation and physical trial history","isoparametric interpolation including old nodal fields","Gauss integration","most residual assembly","complex-step material/H-carrier tangent blocks and most tangent assembly"],
              handwritten_responsibilities=["authored Python physics and mixed formulation","nodal damage-reaction quadrature and exact chain-rule/diagonal assembly extension","dimension/material/geometry/finite checks","transactional state output","DTIME=0 and invalid-state PNEWDT"],
              material_verify_scope='Material/H-carrier derivatives; final nodal reaction is checked separately by independent assembled residual and FD tests.',
              unchanged_host_and_tool_audit=audit,
              element_sha256=element,host_sha256=host,generator_sha256=generator,
              dof_names=contract["dof_names"],history_slots=contract["history_slots"],
              props_order=contract["prop_names"],abaqus_execution=False)
    (ROOT/"verification/provenance.json").write_text(json.dumps(data,indent=2)+"\n")
    print(json.dumps({k:data[k] for k in ("python_version","gfortran_version","fortran_line_allocation","unchanged_host_and_tool_audit")},indent=2))


if __name__=="__main__":
    main()
