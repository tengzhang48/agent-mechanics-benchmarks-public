"""Candidate: accepted-history split with Gauss row-sum reaction lumping."""
from pathlib import Path
import sys
import json

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from verification.exploration.build_lump_gauss import LumpedCoreQuad
from element.model import verification_states
import abaqus_ufl as au


if __name__=="__main__":
    out=ROOT/"tmp/lag_lumped";out.mkdir(exist_ok=True)
    p=LumpedCoreQuad();path=out/"core.for"
    au.generate_uel(p,str(path),element="quad4",formulation="standard",mat_prefix="cf",
        verify=dict(state=verification_states(p)[0],tol=2e-5))
    text=path.read_text().replace("SUBROUTINE UEL(","SUBROUTINE CHEMO_KERNEL(")
    text=text.replace("END SUBROUTINE UEL","END SUBROUTINE CHEMO_KERNEL")
    a=text.index('C       Reject invalid deformation states: det F must stay')
    b=text.index('        END IF',a)+len('        END IF')
    text=text[:a]+'C       Small-strain model: no det(I+grad u) restriction.'+text[b:]
    path.write_text(text)
    wrapper=(ROOT/"verification/exploration/gauss_intact/uel.f90").read_text()
    marker="  ok = all(ieee_is_finite(RHS"
    call="  call accepted_damage(RHS,AMATRX,COORDS,U,saved,PROPS,NDOFEL,MLVARX,PNEWDT)\n"
    (out/"uel.f90").write_text(wrapper.replace(marker,call+marker))
    c=json.loads((ROOT/"verification/exploration/gauss_intact/uel_contract.json").read_text())
    c["source"]=["tmp/lag_lumped/core.for","verification/exploration/accepted_damage.f90","tmp/lag_lumped/uel.f90"]
    (out/"contract.json").write_text(json.dumps(c,indent=2)+"\n")
