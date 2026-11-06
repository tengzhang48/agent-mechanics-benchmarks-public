! Mass-lumped damage reaction, with Gauss-point history from generated cf_trial_history.
! This is the only handwritten FE residual/tangent contribution. The generated
! core supplies mechanics, transport, pressure, alias and damage-gradient terms.
subroutine stable_damage(RHS,K,COORDS,U,OLD,PROPS,NDOFEL,MLVARX,PNEWDT)
  implicit none
  integer,intent(in)::NDOFEL,MLVARX
  double precision,intent(in)::COORDS(2,4),U(NDOFEL),OLD(12),PROPS(11)
  double precision,intent(inout)::RHS(MLVARX,*),K(NDOFEL,NDOFEL),PNEWDT
  double precision::xi(4,2),wg(4),N(4),dn(4,2),B(4,2),Ji(2,2),detj
  double precision::F(3,3),c,s,sa,weight,hh,aa,dh(5),a0
  complex(kind=8)::Fz(3,3),cz,sz,ho,ao,wo,hv(3)
  double precision,parameter::cs=1d-20
  integer::ng,g,a,bn,i,j,col,row,stat
  aa=1d0-PROPS(5)
  call gauss_2d_2x2(xi,wg,ng)
  do g=1,ng
    call shape_quad4(xi(g,1),xi(g,2),N,dn)
    call map_grad_2d(dn,COORDS,4,B,detj,Ji,stat)
    if(stat==0)then
      PNEWDT=0.5d0
      return
    endif
    call eye33d(F)
    c=0d0
    s=0d0
    do a=1,4
      do i=1,2
        do j=1,2
          F(i,j)=F(i,j)+U(6*(a-1)+i)*B(a,j)
        enddo
      enddo
      c=c+U(6*a-3)*N(a)
      s=s+U(6*a)*N(a)
    enddo
    Fz=cmplx(F,0d0,kind=8)
    cz=cmplx(c,0d0,kind=8)
    sz=cmplx(s,0d0,kind=8)
    ho=cmplx(OLD(3*g-2),0d0,kind=8)
    ao=cmplx(OLD(3*g-1),0d0,kind=8)
    wo=cmplx(OLD(3*g),0d0,kind=8)
    call cf_trial_history(Fz,cz,sz,ho,ao,wo,PROPS,hv)
    hh=real(hv(1),kind=8)
    ! Only three strain-gradient components are independent. For this fixed
    ! small-strain model dH/dF12 equals dH/dF21, as checked independently.
    Fz(1,1)=Fz(1,1)+cmplx(0d0,cs,kind=8)
    call cf_trial_history(Fz,cz,sz,ho,ao,wo,PROPS,hv)
    dh(1)=aimag(hv(1))/cs
    Fz=cmplx(F,0d0,kind=8)
    Fz(2,2)=Fz(2,2)+cmplx(0d0,cs,kind=8)
    call cf_trial_history(Fz,cz,sz,ho,ao,wo,PROPS,hv)
    dh(2)=aimag(hv(1))/cs
    Fz=cmplx(F,0d0,kind=8)
    Fz(1,2)=Fz(1,2)+cmplx(0d0,cs,kind=8)
    call cf_trial_history(Fz,cz,sz,ho,ao,wo,PROPS,hv)
    dh(3)=aimag(hv(1))/cs
    Fz=cmplx(F,0d0,kind=8)
    cz=cmplx(c,cs,kind=8)
    call cf_trial_history(Fz,cz,sz,ho,ao,wo,PROPS,hv)
    dh(4)=aimag(hv(1))/cs
    cz=cmplx(c,0d0,kind=8)
    sz=cmplx(s,cs,kind=8)
    call cf_trial_history(Fz,cz,sz,ho,ao,wo,PROPS,hv)
    dh(5)=aimag(hv(1))/cs
    a0=1d0+aa*hh
    do a=1,4
      row=6*a-2
      weight=wg(g)*detj*N(a)
      sa=U(6*a)
      RHS(row,1)=RHS(row,1)-weight*(a0*sa-1d0)
      K(row,6*a)=K(row,6*a)+weight*a0
      do bn=1,4
        col=6*bn-5
        K(row,col)=K(row,col)+weight*aa*sa*(dh(1)*B(bn,1)+dh(3)*B(bn,2))
        K(row,col+1)=K(row,col+1)+weight*aa*sa*(dh(3)*B(bn,1)+dh(2)*B(bn,2))
        K(row,col+2)=K(row,col+2)+weight*aa*sa*dh(4)*N(bn)
        K(row,col+5)=K(row,col+5)+weight*aa*sa*dh(5)*N(bn)
      enddo
    enddo
  enddo
end subroutine stable_damage
