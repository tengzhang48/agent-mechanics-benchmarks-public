"""Exact nodal left preconditioning: scale each damage row by sqrt(g(s_a)).

The factor is identical for a shared node in all neighboring elements, so the
assembled nonlinear roots are unchanged. The tangent includes its derivative.
"""
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[2]

ROW_SCALE="""  ! Nodal left scaling preserves assembled roots because this positive
  ! factor depends only on the shared nodal intactness. Include d(factor)/dU.
  do i=1,4
    row=6*i-2
    col=6*i
    raw_residual=-RHS(row,1)
    factor=sqrt((1d0-PROPS(5))*U(col)**2+PROPS(5))
    dfactor=(1d0-PROPS(5))*U(col)/factor
    AMATRX(row,:)=factor*AMATRX(row,:)
    AMATRX(row,col)=AMATRX(row,col)+dfactor*raw_residual
    RHS(row,1)=factor*RHS(row,1)
  enddo
"""


if __name__=="__main__":
    folder=ROOT/"tmp/row_scaling";folder.mkdir(exist_ok=True)
    wrapper=(ROOT/"verification/exploration/gauss_intact/uel.f90").read_text()
    wrapper=wrapper.replace("  integer :: i","  integer :: i,row,col\n  double precision :: raw_residual,factor,dfactor")
    marker="  ok = all(ieee_is_finite(RHS"
    wrapper=wrapper.replace(marker,ROW_SCALE+marker)
    (folder/"uel.f90").write_text(wrapper)
    c=json.loads((ROOT/"verification/exploration/gauss_intact/uel_contract.json").read_text())
    c["source"]=["verification/exploration/gauss_intact/generated_kernel.for","tmp/row_scaling/uel.f90"]
    (folder/"contract.json").write_text(json.dumps(c,indent=2)+"\n")
