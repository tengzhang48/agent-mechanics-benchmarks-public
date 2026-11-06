"""Known-broken generated controls: verify() can pass incorrect physical laws."""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from element.model import ChemoFatigue,ChemoFatigueQuad,verification_states
from verification.common import compiled
from verification.reference import residual
import fehost
import abaqus_ufl as au


class WrongPressure(ChemoFatigue):
    def pressure_resid(self,F,c,intact,p):
        mu_lame=self.E/(2.0*(1.0+self.nu))
        lam=self.E*self.nu/((1.0+self.nu)*(1.0-2.0*self.nu))
        tr=F[0,0]+F[1,1]-2.0-2.0*self.beta*(c-self.c0)
        degr=(1.0-self.kappa)*intact**2+self.kappa
        return p-degr*2.0*(mu_lame+lam)*tr/2.0  # deliberate /2 defect


class WrongDeviator(ChemoFatigue):
    def _psi_plus(self,F,c):
        mu_lame=self.E/(2.0*(1.0+self.nu))
        lam=self.E*self.nu/((1.0+self.nu)*(1.0-2.0*self.nu))
        bulk=lam+2.0*mu_lame/3.0
        xx=F[0,0]-1.0-self.beta*(c-self.c0)
        yy=F[1,1]-1.0-self.beta*(c-self.c0)
        xy=0.5*(F[0,1]+F[1,0])
        tr=xx+yy
        pos=tr
        if pos.real<0.0:
            pos=0.0
        # Deliberately use the 3D embedded deviator, rather than dev_2.
        return 0.5*bulk*pos**2+mu_lame*((xx-tr/3.0)**2+(yy-tr/3.0)**2+2.0*xy**2+(tr/3.0)**2)


def main():
    archive=ROOT/"verification/exploration/gauss_intact"
    base=json.loads((archive/"uel_contract.json").read_text())
    prm=dict(zip(base["prop_names"],base["props"]))
    X=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
    U=np.zeros((4,6));U[:,:2]=.002*X;U[:,2]=prm["c0"];U[:,3]=.2;U[:,5]=.8
    mu=prm["E"]/(2*(1+prm["nu"]));lam=prm["E"]*prm["nu"]/((1+prm["nu"])*(1-2*prm["nu"]))
    g=(1-prm["kappa"])*.8**2+prm["kappa"]
    U[:,4]=g*4*(mu+lam)*.002/3
    old=np.zeros_like(U);old[:,2]=prm["c0"];old[:,5]=1.
    R0,S0=residual(X,U.ravel(),old.ravel(),np.zeros(12),.2,prm)
    records=[]
    for cls in (WrongPressure,WrongDeviator):
        class MutantProblem(ChemoFatigueQuad):
            material=cls
        problem=MutantProblem()
        assert problem.verify(state=verification_states(problem)[0],verbose=False,tol=2e-5)
        folder=ROOT/"tmp/mutations"/cls.__name__;folder.mkdir(parents=True,exist_ok=True)
        path=folder/"kernel.for"
        au.generate_uel(problem,str(path),element="quad4",formulation="standard",mat_prefix="cf",
            verify=dict(state=verification_states(problem)[0],tol=2e-5))
        txt=path.read_text().replace("SUBROUTINE UEL(","SUBROUTINE CHEMO_KERNEL(")
        path.write_text(txt.replace("END SUBROUTINE UEL","END SUBROUTINE CHEMO_KERNEL"))
        contract={**base,"source":[str(path.relative_to(ROOT)),str((archive/"uel.f90").relative_to(ROOT))]}
        m=compiled(contract)
        el=fehost.ContractUEL(m,contract,1,.2)
        rr,kk,ss,pn=el._call(0,X,U.ravel(),U.ravel()-old.ravel(),np.zeros(12))
        rerr=float(np.max(abs(-rr-R0)));serr=float(np.max(abs(ss-S0)))
        assert rerr>1e-5 or serr>1e-5,(cls.__name__,rerr,serr)
        records.append(dict(mutation=cls.__name__,material_verify_passes=True,
            compiled_independent_oracle_rejects=True,residual_max_error=rerr,state_max_error=serr))
    result=dict(status="passed",controls=records)
    (ROOT/"verification/mutation_checks.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
