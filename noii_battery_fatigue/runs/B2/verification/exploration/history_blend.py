"""Generalized history tangent in a narrow switch band; exact real histories."""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from element.model import ChemoFatigue,ChemoFatigueQuad,verification_states
from abaqus_ufl.core.tensor import log
import abaqus_ufl as au


class HistoryBlend(ChemoFatigue):
    def _trial_history(self,F,c,intact,hist_old,acc_old,wprev_old):
        psi_plus=self._psi_plus(F,c)
        degr=(1.0-self.kappa)*intact**2+self.kappa
        wnow=degr*psi_plus
        anew=acc_old
        if wnow.real>wprev_old.real:
            anew=acc_old+wnow-wprev_old
        fatigue=1.0
        if anew.real>self.alpha_T:
            inner=1.0-self.k_fg*log(anew/self.alpha_T)/log(10.0)
            if inner.real<0.0:
                inner=0.0
            fatigue=inner*inner
            if fatigue.real<1e-15:
                fatigue=1e-15
        candidate=2.0*self.ell*psi_plus/(self.G_c*fatigue)
        hnew=candidate
        if hnew.real<hist_old.real:
            hnew=hist_old
        gap=candidate.real-hist_old.real
        band=1e-8*(1.0+hist_old.real)
        if gap.real<band.real:
            if gap.real>-band.real:
                hnew=hnew.real+0.5*(candidate-candidate.real)
        return np.array([hnew,anew,wnow])


class BlendQuad(ChemoFatigueQuad):
    material=HistoryBlend


if __name__=="__main__":
    folder=ROOT/"tmp/history_blend";folder.mkdir(exist_ok=True)
    problem=BlendQuad()
    path=folder/"core.for"
    au.generate_uel(problem,str(path),element="quad4",formulation="standard",mat_prefix="cf",
        verify=dict(state=verification_states(problem)[0],tol=2e-5))
    txt=path.read_text().replace("SUBROUTINE UEL(","SUBROUTINE CHEMO_KERNEL(")
    path.write_text(txt.replace("END SUBROUTINE UEL","END SUBROUTINE CHEMO_KERNEL"))
    c=json.loads((ROOT/"verification/exploration/gauss_intact/uel_contract.json").read_text())
    c["source"]=["tmp/history_blend/core.for","verification/exploration/gauss_intact/uel.f90"]
    (folder/"contract.json").write_text(json.dumps(c,indent=2)+"\n")
