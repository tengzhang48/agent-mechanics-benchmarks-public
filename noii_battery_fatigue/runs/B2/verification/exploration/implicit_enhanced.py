"""Late controlled alternative: current-H Gauss damage with enhanced strain.

This never changes element/. It tests whether the strain enhancement alone
removes the earlier implicit-history failure. The condensed damage coupling
is obtained by complex-step differentiation of the generated point law.
"""
from pathlib import Path
import sys
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import abaqus_ufl as au
from element.model import ChemoFatigueQuad,verification_states
from verification.exploration.build_enhanced import inject
from verification.reference import residual,shape,GP
from verification.common import compiled,properties
import fehost


def build():
    out=ROOT/'tmp/implicit_enhanced';out.mkdir(exist_ok=True)
    problem=ChemoFatigueQuad()
    path=out/'core.for'
    au.generate_uel(problem,str(path),element='quad4',formulation='standard',mat_prefix='cf',
                    verify=dict(state=verification_states(problem)[0],tol=2e-5))
    text=path.read_text().replace('SUBROUTINE UEL(','SUBROUTINE CHEMO_KERNEL(')
    text=text.replace('END SUBROUTINE UEL','END SUBROUTINE CHEMO_KERNEL')
    start=text.index('C       Reject invalid deformation states: det F must stay')
    end=text.index('        END IF',start)+len('        END IF')
    text=text[:start]+'C       Small-strain model: no det(I+grad u) restriction.'+text[end:]
    text=inject(text).replace('CALL enhanced_modes(coords_2d,U,PROPS,fcorr,kcorr,PNEWDT)',
                             'CALL enhanced_modes(coords_2d,U,PROPS,SVARS,fcorr,kcorr,PNEWDT)')
    path.write_text(text)
    modes=(ROOT/'element/enhanced_modes.f90').read_text()
    modes=modes.replace('subroutine enhanced_modes(coords,U,props,fcorr,kcorr,pnewdt)',
                        'subroutine enhanced_modes(coords,U,props,old,fcorr,kcorr,pnewdt)')
    modes=modes.replace('double precision,intent(in)::coords(2,4),U(24),props(11)',
                        'double precision,intent(in)::coords(2,4),U(24),props(11),old(12)')
    modes=modes.replace('  integer::ng,k,a,bn,i,j,l,col,row,stat,ip,ids(8)',
        '  double complex::fbase(3,3),fz(3,3),phasez,cz,sz,hz,az,wz\n'
        '  double precision::dqphase(4)\n'
        '  integer::ng,k,a,bn,i,j,l,col,row,stat,ip,ids(8)')
    marker='    gp=2d0*(1d0-props(5))*sval(k)'
    code='''
    cg=dot_product(Nall(:,k),U(3:24:6))
    cz=dcmplx(cg,0d0);sz=dcmplx(sval(k),0d0)
    hz=dcmplx(old(3*k-2),0d0);az=dcmplx(old(3*k-1),0d0);wz=dcmplx(old(3*k),0d0)
    fbase=(0d0,0d0)
    fbase(1,1)=1d0+eps(1,k)+h(1)+beta*(cg-props(4))
    fbase(2,2)=1d0+eps(2,k)+h(2)+beta*(cg-props(4))
    fbase(3,3)=1d0
    fbase(1,2)=0.5d0*(eps(3,k)+h(3));fbase(2,1)=fbase(1,2)
    do l=1,4
      fz=fbase
      fz(1,1)=fz(1,1)+(0d0,1d-20)*eq(1,l,k)
      fz(2,2)=fz(2,2)+(0d0,1d-20)*eq(2,l,k)
      fz(1,2)=fz(1,2)+(0d0,5d-21)*eq(3,l,k)
      fz(2,1)=fz(2,1)+(0d0,5d-21)*eq(3,l,k)
      call cf_phase_storage(fz,cz,sz,hz,az,wz,props,phasez)
      dqphase(l)=aimag(phasez)/1d-20
    enddo
'''
    assert modes.count(marker)==1;modes=modes.replace(marker,code+marker)
    marker='      row=6*a-1'
    modes=modes.replace(marker,'      xq(6*a-2,:)=xq(6*a-2,:)+wt(k)*Nall(a,k)*dqphase\n'+marker)
    (out/'modes.f90').write_text(modes)
    guard=(ROOT/'element/uel.f90').read_text()
    line='  call accepted_damage(RHS,AMATRX,COORDS,U,saved,PROPS,NDOFEL,MLVARX,PNEWDT)\n'
    assert line in guard;guard=guard.replace(line,'')
    (out/'guard.f90').write_text(guard)
    c=json.loads((ROOT/'element/uel_contract.json').read_text())
    c.update(source=[str(p.relative_to(ROOT)) for p in (path,out/'modes.f90',out/'guard.f90')],
        time_scheme='Backward Euler and current endpoint H',damage_reaction='Gauss current H',
        formulation_id='implicit_enhanced_experiment')
    (out/'contract.json').write_text(json.dumps(c,indent=2)+'\n')
    return c


def oracle(X,U,prev,old,dt,p):
    r,state=residual(X,U,prev,old,dt,p,rule='enhanced_lagged')
    v=U.reshape(4,6);q=old.reshape(4,3);a=1-p['kappa']
    for k,(xi,eta) in enumerate(GP):
        N,dn=shape(xi,eta);det=np.linalg.det(X.T@dn);sg=N@v[:,5]
        r[3::6]+=N*((1+a*state[3*k])*sg-(1+a*q[k,0])*v[:,5])*det
    return r,state


def check(c):
    el=fehost.ContractUEL(compiled(c,checked=True),c,1,.2);p=properties(c)
    X=np.array([[0.,0.],[1.2,.1],[1.03,.91],[-.08,1.05]])
    prev=np.zeros((4,6));prev[:,2]=p['c0'];prev[:,5]=1.
    v=prev.copy();v[:,:2]=X@np.array([[.005,.001],[.002,-.001]]).T
    v[:,2]+=[.02,.07,.10,.04];v[:,3]=[.13,.35,.48,.26];v[:,5]=1-v[:,3];v[:,4]=[.02,.03,-.02,-.01]
    U=v.ravel();pr=prev.ravel();old=np.array([[.5,2e-4,1e-6],[.2,4e-4,2e-6],[20.,5e-4,.003],[2.,7e-4,1e-6]]).ravel()
    rhs,K,s,pn=el._call(0,X,U,U-pr,old.copy())
    R=-rhs;rr,ss=oracle(X,U,pr,old,.2,p)
    assert pn==1 and np.allclose(R,rr,rtol=2e-10,atol=2e-11)
    assert np.allclose(s,ss,rtol=2e-10,atol=2e-11)
    fd=np.empty_like(K);eps=3e-8
    for j in range(24):
        step=np.zeros(24);step[j]=eps
        a=el._call(0,X,U+step,U+step-pr,old.copy())[0]
        b=el._call(0,X,U-step,U-step-pr,old.copy())[0]
        fd[:,j]=(-a+b)/(2*eps)
    err=float(np.max(abs(K-fd))/np.max(abs(K)));assert err<2e-6,err
    result=dict(status='passed',residual_max_error=float(np.max(abs(R-rr))),
        state_max_error=float(np.max(abs(s-ss))),tangent_max_normalized_error=err,
        scope='Current-H enhanced Gauss alternative, independent non-affine residual and full condensed tangent')
    (ROOT/'verification/exploration/implicit_enhanced_checks.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':check(build())
