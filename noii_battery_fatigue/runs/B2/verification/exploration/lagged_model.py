"""First-order history splitting candidate (real endpoint history kept exact).

Damage uses accepted H_n; history H_{n+1},alpha_{n+1},W_{n+1} is recomputed
from the converged fields for the next step. This is a time discretization
with an explicit history lag, and must be assessed under time refinement.
"""
from pathlib import Path
import sys
import json

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from element.model import ChemoFatigue,ChemoFatigueQuad,verification_states
import abaqus_ufl as au


class LaggedMaterial(ChemoFatigue):
    def phase_storage(self,intact,hist_old):
        return (1.0+(1.0-self.kappa)*hist_old)*intact-1.0


class LaggedQuad(ChemoFatigueQuad):
    material=LaggedMaterial

    def phase_equation(self,eta,intact,grad_intact,hist_old):
        return (self.material.phase_storage(intact,hist_old),self.material.phase_flux(grad_intact))


if __name__=="__main__":
    out=ROOT/"tmp/lagged";out.mkdir(exist_ok=True)
    p=LaggedQuad()
    path=out/"core.for"
    au.generate_uel(p,str(path),element="quad4",formulation="standard",mat_prefix="cf",
        verify=dict(state=verification_states(p)[0],tol=2e-5))
    text=path.read_text().replace("SUBROUTINE UEL(","SUBROUTINE CHEMO_KERNEL(")
    path.write_text(text.replace("END SUBROUTINE UEL","END SUBROUTINE CHEMO_KERNEL"))
    c=json.loads((ROOT/"verification/exploration/gauss_intact/uel_contract.json").read_text())
    c["source"]=["tmp/lagged/core.for","verification/exploration/gauss_intact/uel.f90"]
    (out/"contract.json").write_text(json.dumps(c,indent=2)+"\n")
