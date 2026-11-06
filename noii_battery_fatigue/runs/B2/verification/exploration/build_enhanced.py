"""Experimental extension of the selected generated kernel by four local modes."""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[2]


def inject(text):
    declaration="      DOUBLE PRECISION :: fcorr(3,3,4), kcorr(24,24)\n"
    marker="      INTEGER :: idx, ii_v, jj_v, i, j, k, l, kk, row, col"
    assert text.count(marker)==1
    text=text.replace(marker,declaration+marker)
    marker="C     Gauss quadrature\n"
    code="      CALL enhanced_modes(coords_2d,U,PROPS,fcorr,kcorr,PNEWDT)\n      IF (PNEWDT .LT. 1.0D0) RETURN\n\n"
    text=text.replace(marker,code+marker,1)
    marker="C       Small-strain model: no det(I+grad u) restriction."
    assert text.count(marker)==1
    text=text.replace(marker,"        F = F + fcorr(:,:,kk)\n"+marker)
    marker="C     --- Non-finite RHS or AMATRX: request a smaller increment"
    text=text.replace(marker,"      AMATRX = AMATRX - kcorr\n\n"+marker,1)
    return text


if __name__=="__main__":
    out=ROOT/"tmp/enhanced";out.mkdir(exist_ok=True)
    # This predecessor was promoted into element/. Reproduce the experiment
    # from its preserved, unenriched input rather than injecting it twice.
    text=inject((ROOT/"verification/exploration/basic_split_corrected/generated_kernel.for").read_text())
    (out/"core.for").write_text(text)
    c=json.loads((ROOT/"element/uel_contract.json").read_text())
    c["source"]=["tmp/enhanced/core.for","verification/exploration/enhanced_modes.f90",
                 "element/accepted_damage.f90","element/uel.f90"]
    c["element_local_variables"]=["ux_xi2","ux_eta2","uy_xi2","uy_eta2"]
    c["local_variable_treatment"]="Exact static condensation of a 4x4 positive-definite elastic system; recomputed every call"
    (out/"contract.json").write_text(json.dumps(c,indent=2)+"\n")
