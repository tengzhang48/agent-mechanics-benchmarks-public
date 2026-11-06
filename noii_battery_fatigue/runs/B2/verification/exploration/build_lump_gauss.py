"""Candidate: keep Gauss history, lump only the damage reaction matrix."""
from pathlib import Path
import sys
import json

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from element.model import ChemoFatigue,ChemoFatigueQuad,verification_states
import abaqus_ufl as au


class LumpedCoreMaterial(ChemoFatigue):
    def phase_storage(self):
        return 0.0


class LumpedCoreQuad(ChemoFatigueQuad):
    material=LumpedCoreMaterial

    def phase_equation(self,eta,grad_intact):
        return self.material.phase_storage(),self.material.phase_flux(grad_intact)


def generate():
    folder=ROOT/"tmp/lump_gauss";folder.mkdir(exist_ok=True)
    p=LumpedCoreQuad()
    for s in verification_states(p):
        p.verify(state=s,verbose=False,tol=2e-5)
    path=folder/"core.for"
    au.generate_uel(p,str(path),element="quad4",formulation="standard",mat_prefix="cf",
        verify=dict(state=verification_states(p)[0],tol=2e-5))
    text=path.read_text().replace("SUBROUTINE UEL(","SUBROUTINE CHEMO_KERNEL(")
    path.write_text(text.replace("END SUBROUTINE UEL","END SUBROUTINE CHEMO_KERNEL"))
    wrapper=(ROOT/"verification/exploration/gauss_intact/uel.f90").read_text()
    marker="  ok = all(ieee_is_finite(RHS"
    call="  call stable_damage(RHS,AMATRX,COORDS,U,saved,PROPS,NDOFEL,MLVARX,PNEWDT)\n"
    wrapper=wrapper.replace(marker,call+marker)
    (folder/"uel.f90").write_text(wrapper)
    c=json.loads((ROOT/"verification/exploration/gauss_intact/uel_contract.json").read_text())
    c["source"]=["tmp/lump_gauss/core.for","verification/exploration/lump_gauss_damage.f90","tmp/lump_gauss/uel.f90"]
    (folder/"contract.json").write_text(json.dumps(c,indent=2)+"\n")


if __name__=="__main__":
    generate()
