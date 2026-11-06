"""Regenerate the same Q4 mixed formulation with 3x3 quadrature for comparison.

Uses the public ElementConfig API; it does not modify the generator or host.
The alternate source, contract and compiler products all stay under tmp/.
"""
from dataclasses import replace
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

from common import ROOT, fehost
from run_case import run

sys.path[:0]=[str(ROOT/"tools/abaqus_ufl"),str(ROOT/"element")]
from abaqus_ufl import ELEMENT_CONFIGS, generate_uel
from model import ChemoProblem, verification_states


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--version',default='selected')
    args=parser.parse_args()
    alt=ROOT/"tmp/quadrature9"/args.version
    (alt/"element").mkdir(parents=True,exist_ok=True)
    problem=ChemoProblem()
    cfg=replace(ELEMENT_CONFIGS["quad4"],n_gauss_points=9,gauss_subroutine="gauss_2d_3x3")
    path=alt/"element/generated_model.for"
    generate_uel(problem,str(path),element_config=cfg,formulation="standard",mat_prefix="chem",
                 verify=dict(state=next(verification_states())[1]))
    text=path.read_text().replace("SUBROUTINE UEL(","SUBROUTINE CHEM_UEL_CORE(",1)
    text=text.replace("END SUBROUTINE UEL","END SUBROUTINE CHEM_UEL_CORE")
    if (ROOT/'element/lump_reaction.py').exists():
        from lump_reaction import transform
        text,_=transform(text)
    path.write_text(text)
    contract=fehost.load_contract(str(ROOT))
    old_nsv=contract['nsvars']
    stride=old_nsv//len(contract['history_slots']['H'])
    new_nsv=9*stride
    wrapper=(ROOT/"element/uel.f90").read_text()
    wrapper=wrapper.replace(f"trial_state({old_nsv})",f"trial_state({new_nsv})")
    wrapper=wrapper.replace(f"NSVARS /= {old_nsv}",f"NSVARS /= {new_nsv}")
    wrapper=wrapper.replace(f"do i = 1, {old_nsv}, {stride}",f"do i = 1, {new_nsv}, {stride}")
    (alt/"element/uel.f90").write_text(wrapper)
    contract["nsvars"]=new_nsv
    contract["history_slots"]={"H":list(range(1,new_nsv+1,stride)),"alpha":list(range(2,new_nsv+1,stride))}
    (alt/"element/uel_contract.json").write_text(json.dumps(contract,indent=2)+"\n")
    manifest={name:hashlib.sha256((alt/name).read_bytes()).hexdigest()
              for name in contract['source']+['element/uel_contract.json']}
    snapshot=ROOT/f'verification/quadrature_comparison/{args.version}/source_snapshot'
    (snapshot/'element').mkdir(parents=True,exist_ok=True)
    for name in manifest:
        shutil.copyfile(alt/name,snapshot/name)
    (snapshot/'provenance.json').write_text(json.dumps(dict(
        n_gauss_points=9,element_source_sha256=manifest,
        generation_script='verification/quadrature_study.py',
        base_model_sha256=hashlib.sha256((ROOT/'element/model.py').read_bytes()).hexdigest(),
        base_extension_sha256=hashlib.sha256((ROOT/'element/lump_reaction.py').read_bytes()).hexdigest()
        if (ROOT/'element/lump_reaction.py').exists() else None),indent=2)+'\n')
    tag=hashlib.sha256(json.dumps(manifest,sort_keys=True).encode()).hexdigest()[:12]
    module=fehost.build_module(str(alt),contract,str(alt/"tmp/build"),"uel_quad9_"+tag)
    for n in (16,32):
        out=ROOT/f"verification/quadrature_comparison/{args.version}/n{n}_3x3"
        summary=run(ROOT/"host/dev_problem/d0.json",n=n,cycles=3,spc=80,out=out,
                    module_contract=(module,contract),source_manifest=manifest)
        assert summary["status"]=="ok"


if __name__=="__main__":
    main()
