! First-order accepted-history split, with a row-sum lumped reaction matrix.
! H is from the START of this increment. Generated history updates are returned
! to the host separately and become available only after acceptance.
subroutine accepted_damage(RHS,K,COORDS,U,OLD,PROPS,NDOFEL,MLVARX,PNEWDT)
  implicit none
  integer,intent(in)::NDOFEL,MLVARX
  double precision,intent(in)::COORDS(2,4),U(NDOFEL),OLD(12),PROPS(11)
  double precision,intent(inout)::RHS(MLVARX,*),K(NDOFEL,NDOFEL),PNEWDT
  double precision::xi(4,2),wg(4),N(4),dn(4,2),B(4,2),Ji(2,2),detj,weight,coef
  integer::ng,g,a,row,stat
  call gauss_2d_2x2(xi,wg,ng)
  do g=1,ng
    call shape_quad4(xi(g,1),xi(g,2),N,dn)
    call map_grad_2d(dn,COORDS,4,B,detj,Ji,stat)
    if(stat==0)then
      PNEWDT=0.5d0
      return
    endif
    coef=1d0+(1d0-PROPS(5))*OLD(3*g-2)
    do a=1,4
      row=6*a-2
      weight=wg(g)*detj*N(a)
      RHS(row,1)=RHS(row,1)-weight*(coef*U(6*a)-1d0)
      K(row,6*a)=K(row,6*a)+weight*coef
    enddo
  enddo
end subroutine accepted_damage
