"""Alternative: half exact bulk stress, half projected hydrostatic stress.

The compliance-weighted pressure constraint keeps the continuum constitutive
law unchanged. A nonzero exact bulk fraction keeps the local mode system SPD.
"""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from element.model import LumpedCoreMaterial,ChemoFatigueQuad,verification_states
from verification.exploration.build_enhanced import inject
import abaqus_ufl as au


class HybridMaterial(LumpedCoreMaterial):
    def _make_state(self,**kwargs):
        kwargs.setdefault("p",0.)
        return super()._make_state(**kwargs)

    def stress_PK1(self,F,c,intact,p,hist_old,acc_old,wprev_old,dt):
        mu_lame=self.E/(2.0*(1.0+self.nu))
        lam=self.E*self.nu/((1.0+self.nu)*(1.0-2.0*self.nu))
        xx=F[0,0]-1.0-self.beta*(c-self.c0)
        yy=F[1,1]-1.0-self.beta*(c-self.c0)
        xy=0.5*(F[0,1]+F[1,0])
        g=(1.0-self.kappa)*intact**2+self.kappa
        mean=0.5*g*(mu_lame+lam)*(xx+yy)+0.75*p
        sxx=g*mu_lame*(xx-yy)+mean
        syy=-g*mu_lame*(xx-yy)+mean
        sxy=g*2.0*mu_lame*xy
        stress=np.array([[sxx,sxy,0.0],[sxy,syy,0.0],[0.0,0.0,0.0]])
        hv=self._trial_history(F,c,intact,hist_old,acc_old,wprev_old)
        return stress,{"hist":hv[0],"acc":hv[1],"wprev":hv[2]}

    def pressure_resid(self,F,c,intact,p):
        mu_lame=self.E/(2.0*(1.0+self.nu))
        lam=self.E*self.nu/((1.0+self.nu)*(1.0-2.0*self.nu))
        cp=2.0*(mu_lame+lam)/3.0
        g=(1.0-self.kappa)*intact**2+self.kappa
        tr=F[0,0]+F[1,1]-2.0-2.0*self.beta*(c-self.c0)
        return p/(g*cp)-tr


class HybridQuad(ChemoFatigueQuad):
    material=HybridMaterial

    def momentum_equation(self,v,F,c,intact,p,hist_old,acc_old,wprev_old,dt):
        return self.material.stress_PK1(F,c,intact,p,hist_old,acc_old,wprev_old,dt)

    def phase_equation(self,eta,grad_intact):
        return self.material.phase_storage(),self.material.phase_flux(grad_intact)


if __name__=="__main__":
    out=ROOT/"tmp/hybrid_enhanced";out.mkdir(exist_ok=True)
    p=HybridQuad();path=out/"core.for"
    au.generate_uel(p,str(path),element="quad4",formulation="standard",mat_prefix="cf",
        verify=dict(state=verification_states(p)[0],tol=2e-5))
    text=path.read_text().replace("SUBROUTINE UEL(","SUBROUTINE CHEMO_KERNEL(")
    text=text.replace("END SUBROUTINE UEL","END SUBROUTINE CHEMO_KERNEL")
    a=text.index('C       Reject invalid deformation states: det F must stay')
    b=text.index('        END IF',a)+len('        END IF')
    text=text[:a]+'C       Small-strain model: no det(I+grad u) restriction.'+text[b:]
    path.write_text(inject(text))
    mode=(ROOT/"verification/exploration/enhanced_modes.f90").read_text()
    mode=mode.replace('double precision::eps(3,4)', 'double precision::pv(4),bulk,tau\n  double precision::eps(3,4)')
    mode=mode.replace('CE=0d0;CE(1,1)=lam+2d0*mu;CE(2,2)=CE(1,1)\n  CE(1,2)=lam;CE(2,1)=lam;CE(3,3)=mu',
        'tau=0.5d0;bulk=lam+mu\n  CE=0d0;CE(1,1)=mu+tau*bulk;CE(2,2)=CE(1,1)\n  CE(1,2)=-mu+tau*bulk;CE(2,1)=CE(1,2);CE(3,3)=mu')
    mode=mode.replace('cg=0d0;sval(k)=0d0','cg=0d0;sval(k)=0d0;pv(k)=0d0')
    mode=mode.replace('sval(k)=sval(k)+Nall(a,k)*U(6*a)','sval(k)=sval(k)+Nall(a,k)*U(6*a)\n      pv(k)=pv(k)+Nall(a,k)*U(6*a-1)')
    mode=mode.replace('gval(k)*matmul(CE,eps(:,k)))','gval(k)*matmul(CE,eps(:,k))+1.5d0*(1d0-tau)*pv(k)*[1d0,1d0,0d0])')
    mode=mode.replace('Nall(a,k)*gval(k)*cp*(eq(1,:,k)+eq(2,:,k))','Nall(a,k)*(eq(1,:,k)+eq(2,:,k))')
    mode=mode.replace('col=6*a-3\n', 'col=6*a-1\n      qx(:,col)=qx(:,col)+wt(k)*Nall(a,k)*1.5d0*(1d0-tau)*(eq(1,:,k)+eq(2,:,k))\n      col=6*a-3\n')
    (out/"modes.f90").write_text(mode)
    c=json.loads((ROOT/"element/uel_contract.json").read_text())
    c["source"]=["tmp/hybrid_enhanced/core.for","tmp/hybrid_enhanced/modes.f90","element/accepted_damage.f90","element/uel.f90"]
    c["element_local_variables"]=["ux_xi2","ux_eta2","uy_xi2","uy_eta2"]
    c["local_variable_treatment"]="Four zero-mean enhanced strain modes; exact static condensation, recomputed each call"
    (out/"contract.json").write_text(json.dumps(c,indent=2)+"\n")
