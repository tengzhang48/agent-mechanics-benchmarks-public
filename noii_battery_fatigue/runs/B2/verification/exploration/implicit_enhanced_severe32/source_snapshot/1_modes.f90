! Four element-local incompatible displacement modes, statically condensed.
! Their physical gradients have zero volume average, preserving affine patches
! on non-affine Q4 geometry. No local variable is retained between calls.
subroutine enhanced_modes(coords,U,props,old,fcorr,kcorr,pnewdt)
  implicit none
  double precision,intent(in)::coords(2,4),U(24),props(11),old(12)
  double precision,intent(out)::fcorr(3,3,4),kcorr(24,24)
  double precision,intent(inout)::pnewdt
  double precision::xi(4,2),wg(4),Nall(4,4),Ball(4,2,4),wt(4),gm(2,2,4),meanG(2,2)
  double precision::N(4),dn(4,2),B(4,2),Ji(2,2),detj,vol,CE(3,3),eq(3,4,4),bu(3,8,4)
  double precision::eps(3,4),sigma(3),gval(4),sval(4),alocal(4,4),ainv(4,4),rq(4),q(4)
  double precision::xq(24,4),qx(4,24),D(3,4),DC(3),mu,lam,cp,beta,sw,cg,tmp4(4),tmp8(8)
  double precision::aug(4,8),scale,pivot,temp(8),h(3),gp
  double complex::fbase(3,3),fz(3,3),phasez,cz,sz,hz,az,wz
  double precision::dqphase(4)
  integer::ng,k,a,bn,i,j,l,col,row,stat,ip,ids(8)
  fcorr=0d0;kcorr=0d0
  mu=props(1)/(2d0*(1d0+props(2)))
  lam=props(1)*props(2)/((1d0+props(2))*(1d0-2d0*props(2)))
  beta=props(3);cp=2d0*(mu+lam)/3d0
  CE=0d0;CE(1,1)=lam+2d0*mu;CE(2,2)=CE(1,1)
  CE(1,2)=lam;CE(2,1)=lam;CE(3,3)=mu
  DC=matmul(CE,[-beta,-beta,0d0])
  call gauss_2d_2x2(xi,wg,ng)
  meanG=0d0;vol=0d0
  do k=1,4
    call shape_quad4(xi(k,1),xi(k,2),N,dn)
    call map_grad_2d(dn,coords,4,B,detj,Ji,stat)
    if(stat==0)then
      pnewdt=0.5d0
      return
    endif
    wt(k)=wg(k)*detj;Nall(:,k)=N;Ball(:,:,k)=B
    gm(:,1,k)=-2d0*xi(k,1)*Ji(:,1)
    gm(:,2,k)=-2d0*xi(k,2)*Ji(:,2)
    meanG=meanG+wt(k)*gm(:,:,k);vol=vol+wt(k)
  enddo
  meanG=meanG/vol
  eq=0d0;bu=0d0;eps=0d0
  alocal=0d0;rq=0d0
  do k=1,4
    gm(:,:,k)=gm(:,:,k)-meanG
    do i=1,2
      eq(:,i,k)=[gm(1,i,k),0d0,gm(2,i,k)]
      eq(:,i+2,k)=[0d0,gm(2,i,k),gm(1,i,k)]
    enddo
    cg=0d0;sval(k)=0d0
    do a=1,4
      bu(:,2*a-1,k)=[Ball(a,1,k),0d0,Ball(a,2,k)]
      bu(:,2*a,k)=[0d0,Ball(a,2,k),Ball(a,1,k)]
      eps(:,k)=eps(:,k)+bu(:,2*a-1,k)*U(6*a-5)+bu(:,2*a,k)*U(6*a-4)
      cg=cg+Nall(a,k)*U(6*a-3)
      sval(k)=sval(k)+Nall(a,k)*U(6*a)
    enddo
    sw=beta*(cg-props(4));eps(1:2,k)=eps(1:2,k)-sw
    gval(k)=(1d0-props(5))*sval(k)**2+props(5)
    D=gval(k)*matmul(CE,eq(:,:,k))
    alocal=alocal+wt(k)*matmul(transpose(eq(:,:,k)),D)
    rq=rq+wt(k)*matmul(transpose(eq(:,:,k)),gval(k)*matmul(CE,eps(:,k)))
  enddo
  ! Small pivoted inverse; A is the positive-definite local elastic matrix.
  aug=0d0;aug(:,1:4)=alocal
  do i=1,4
    aug(i,4+i)=1d0
  enddo
  scale=maxval(abs(alocal))
  do i=1,4
    ip=i
    do j=i+1,4
      if(abs(aug(j,i))>abs(aug(ip,i)))ip=j
    enddo
    pivot=aug(ip,i)
    if(abs(pivot)<1d-12*scale)then
      pnewdt=0.5d0
      return
    endif
    temp=aug(i,:);aug(i,:)=aug(ip,:);aug(ip,:)=temp
    aug(i,:)=aug(i,:)/aug(i,i)
    do j=1,4
      if(j/=i)aug(j,:)=aug(j,:)-aug(j,i)*aug(i,:)
    enddo
  enddo
  ainv=aug(:,5:8);q=-matmul(ainv,rq)
  xq=0d0;qx=0d0
  do k=1,4
    h=matmul(eq(:,:,k),q)
    fcorr(1,1,k)=h(1);fcorr(2,2,k)=h(2)
    fcorr(1,2,k)=0.5d0*h(3);fcorr(2,1,k)=0.5d0*h(3)
    sigma=matmul(CE,eps(:,k)+h)
    D=gval(k)*matmul(CE,eq(:,:,k))

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
    gp=2d0*(1d0-props(5))*sval(k)
    do a=1,4
      do i=1,2
        row=6*(a-1)+i
        xq(row,:)=xq(row,:)+wt(k)*matmul(bu(:,2*(a-1)+i,k),D)
        qx(:,row)=qx(:,row)+wt(k)*matmul(transpose(D),bu(:,2*(a-1)+i,k))
      enddo
      xq(6*a-2,:)=xq(6*a-2,:)+wt(k)*Nall(a,k)*dqphase
      row=6*a-1
      xq(row,:)=xq(row,:)-wt(k)*Nall(a,k)*gval(k)*cp*(eq(1,:,k)+eq(2,:,k))
      col=6*a-3
      qx(:,col)=qx(:,col)+wt(k)*Nall(a,k)*gval(k)*matmul(transpose(eq(:,:,k)),DC)
      col=6*a
      qx(:,col)=qx(:,col)+wt(k)*Nall(a,k)*gp*matmul(transpose(eq(:,:,k)),sigma)
    enddo
  enddo
  kcorr=matmul(xq,matmul(ainv,qx))
end subroutine enhanced_modes
