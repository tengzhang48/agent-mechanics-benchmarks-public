"""Alternative considered: energy-compatible mixed hydrostatic stress.

Not selected for the deliverable unless its numerical evidence warrants it.
Generate this candidate with this file; do not alter the package or host.
"""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from element.model import ChemoFatigue,ChemoFatigueQuad,verification_states
import abaqus_ufl as au


class MixedStress(ChemoFatigue):
    def _make_state(self,**kwargs):
        kwargs.setdefault("p",0.)
        return super()._make_state(**kwargs)

    def stress_PK1(self,F,c,intact,p,hist_old,acc_old,wprev_old,dt):
        mu_lame=self.E/(2.0*(1.0+self.nu))
        diff=F[0,0]-F[1,1]
        shear=F[0,1]+F[1,0]
        degr=(1.0-self.kappa)*intact**2+self.kappa
        sxx=degr*mu_lame*diff+1.5*p
        syy=-degr*mu_lame*diff+1.5*p
        sxy=degr*mu_lame*shear
        stress=np.array([[sxx,sxy,0.0],[sxy,syy,0.0],[0.0,0.0,0.0]])
        hv=self._trial_history(F,c,intact,hist_old,acc_old,wprev_old)
        return stress,{"hist":hv[0],"acc":hv[1],"wprev":hv[2]}

    def pressure_resid(self,F,c,intact,p):
        mu_lame=self.E/(2.0*(1.0+self.nu))
        lam=self.E*self.nu/((1.0+self.nu)*(1.0-2.0*self.nu))
        factor=2.0*(mu_lame+lam)/3.0
        degr=(1.0-self.kappa)*intact**2+self.kappa
        tr=F[0,0]+F[1,1]-2.0-2.0*self.beta*(c-self.c0)
        return p/(degr*factor)-tr


class MixedQuad(ChemoFatigueQuad):
    material=MixedStress

    def momentum_equation(self,v,F,c,intact,p,hist_old,acc_old,wprev_old,dt):
        return self.material.stress_PK1(F,c,intact,p,hist_old,acc_old,wprev_old,dt)


if __name__=="__main__":
    problem=MixedQuad()
    for state in verification_states(problem):
        problem.verify(state=state,verbose=False,tol=2e-5)
    out=ROOT/"tmp/mixed_stress"
    out.mkdir(exist_ok=True)
    path=out/"kernel.for"
    au.generate_uel(problem,str(path),element="quad4",formulation="standard",mat_prefix="cf",
                    verify=dict(state=verification_states(problem)[0],tol=2e-5))
    text=path.read_text().replace("SUBROUTINE UEL(","SUBROUTINE CHEMO_KERNEL(")
    text=text.replace("END SUBROUTINE UEL","END SUBROUTINE CHEMO_KERNEL")
    path.write_text(text)
    contract=json.loads((ROOT/"verification/exploration/gauss_intact/uel_contract.json").read_text())
    contract["source"]=["tmp/mixed_stress/kernel.for","verification/exploration/gauss_intact/uel.f90"]
    (out/"contract.json").write_text(json.dumps(contract,indent=2)+"\n")
