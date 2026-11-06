"""Candidate split quadrature: Gauss bulk and vertex history/reaction.

Uses two ordinary generated kernels and a handwritten corner quadrature rule.
The history kernel shares the bulk kernel's identical shipped helper library.
"""
from pathlib import Path
from dataclasses import replace
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from element.model import ChemoFatigue,ChemoFatigueQuad,verification_states
import abaqus_ufl as au


class BulkMaterial(ChemoFatigue):
    state_vars={}

    def stress_PK1(self,F,c,intact):
        mu_lame=self.E/(2.0*(1.0+self.nu))
        lam=self.E*self.nu/((1.0+self.nu)*(1.0-2.0*self.nu))
        xx=F[0,0]-1.0-self.beta*(c-self.c0)
        yy=F[1,1]-1.0-self.beta*(c-self.c0)
        xy=0.5*(F[0,1]+F[1,0])
        g=(1.0-self.kappa)*intact**2+self.kappa
        sxx=g*(2.0*mu_lame*xx+lam*(xx+yy))
        syy=g*(2.0*mu_lame*yy+lam*(xx+yy))
        sxy=g*2.0*mu_lame*xy
        return np.array([[sxx,sxy,0.0],[sxy,syy,0.0],[0.0,0.0,0.0]])

    def phase_storage(self):
        return 0.0


class BulkQuad(ChemoFatigueQuad):
    material=BulkMaterial

    def momentum_equation(self,v,F,c,intact):
        return self.material.stress_PK1(F,c,intact)

    def phase_equation(self,eta,grad_intact):
        return self.material.phase_storage(),self.material.phase_flux(grad_intact)


class HistoryMaterial(ChemoFatigue):
    pressure_resid=None
    solvent_flux=None
    solvent_storage=None
    species_flux=None
    species_storage=None

    def stress_PK1(self,F,c,intact,hist_old,acc_old,wprev_old,dt):
        hv=self._trial_history(F,c,intact,hist_old,acc_old,wprev_old)
        stress=np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]])
        return stress,{"hist":hv[0],"acc":hv[1],"wprev":hv[2]}

    def phase_flux(self,grad_intact):
        return 0.0*grad_intact


class HistoryQuad(au.WeakForm):
    material=HistoryMaterial
    ndim=2
    define_fields=ChemoFatigueQuad.define_fields
    momentum_equation=ChemoFatigueQuad.momentum_equation
    phase_equation=ChemoFatigueQuad.phase_equation


CORNER_RULE="""
subroutine corner_quad4(xi,w,ngp)
  implicit none
  integer,intent(out)::ngp
  double precision,intent(out)::xi(4,2),w(4)
  ngp=4
  xi(:,1)=[-1d0,1d0,1d0,-1d0]
  xi(:,2)=[-1d0,-1d0,1d0,1d0]
  w=1d0
end subroutine corner_quad4
"""


def generate():
    out=ROOT/"tmp/lumped_candidate";out.mkdir(exist_ok=True)
    bulk=BulkQuad();history=HistoryQuad()
    for p in (bulk,history):
        for state in verification_states(p):
            p.verify(state=state,verbose=False,tol=2e-5)
    cfg=replace(au.ELEMENT_CONFIGS["quad4"],gauss_subroutine="corner_quad4")
    for name,problem,config in (("bulk",bulk,au.ELEMENT_CONFIGS["quad4"]),("history",history,cfg)):
        path=out/f"{name}.for"
        au.generate_uel(problem,str(path),element_config=config,formulation="standard",
            mat_prefix="c"+name,verify=dict(state=verification_states(problem)[0],tol=2e-5))
        txt=path.read_text().replace("SUBROUTINE UEL(",f"SUBROUTINE CHEMO_{name.upper()}(")
        txt=txt.replace("END SUBROUTINE UEL",f"END SUBROUTINE CHEMO_{name.upper()}")
        if name=="history":
            txt=txt[:txt.index("C     Section 4:")]
        path.write_text(txt)
    wrapper=(ROOT/"verification/exploration/gauss_intact/uel.f90").read_text()
    wrapper=wrapper.replace("  double precision :: saved(NSVARS)",
        "  double precision :: saved(NSVARS), rb(MLVARX,NRHS), kb(NDOFEL,NDOFEL)")
    oldcall=wrapper[wrapper.index("  call CHEMO_KERNEL("):wrapper.index("  ok = all(ieee_is_finite(RHS")]
    newcall=oldcall.replace("CHEMO_KERNEL","CHEMO_BULK")+"  rb=RHS(1:MLVARX,1:NRHS)\n  kb=AMATRX\n"
    newcall+=oldcall.replace("CHEMO_KERNEL","CHEMO_HISTORY")
    newcall+="  RHS(1:MLVARX,1:NRHS)=RHS(1:MLVARX,1:NRHS)+rb\n  AMATRX=AMATRX+kb\n"
    wrapper=wrapper.replace(oldcall,newcall)+CORNER_RULE
    (out/"uel.f90").write_text(wrapper)
    c=json.loads((ROOT/"verification/exploration/gauss_intact/uel_contract.json").read_text())
    c["source"]=["tmp/lumped_candidate/bulk.for","tmp/lumped_candidate/history.for","tmp/lumped_candidate/uel.f90"]
    (out/"contract.json").write_text(json.dumps(c,indent=2)+"\n")


if __name__=="__main__":
    generate()
