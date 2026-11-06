! Direct implementation of the dimensionless model in PUBLIC_TASK.md.
! No persistent memory: the caller owns the transaction on SVARS.
! Six Q1 fields: ux, uy, c, d, q=g*p_eff, z=(1-kappa)*H*(1-d).
! Three-point Gauss integration in each parent direction.
subroutine uel(rhs, amatrx, svars, energy, ndofel, nrhs, nsvars, &
               props, nprops, coords, mcrd, nnode, u, du, v, acc, &
               jtype, time, dtime, kstep, kinc, jelem, params, ndload, &
               jdltyp, adlmag, predef, npredf, lflags, mlvarx, ddlmag, &
               mdload, pnewdt, jprops, njprop, period)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, intent(in) :: ndofel, nrhs, nsvars, nprops, mcrd, nnode
  integer, intent(in) :: jtype, kstep, kinc, jelem, ndload, npredf
  integer, intent(in) :: mlvarx, mdload, njprop
  integer, intent(in) :: jdltyp(mdload,*), lflags(*), jprops(*)
  double precision, intent(in) :: props(nprops), coords(mcrd,nnode)
  double precision, intent(in) :: u(ndofel), du(mlvarx,*), v(ndofel), acc(ndofel)
  double precision, intent(in) :: time(2), dtime, params(*), adlmag(mdload,*)
  double precision, intent(in) :: predef(2,npredf,nnode), ddlmag(mdload,*), period
  double precision, intent(inout) :: svars(nsvars), pnewdt
  double precision, intent(out) :: rhs(mlvarx,nrhs), amatrx(ndofel,ndofel), energy(8)
  integer, parameter :: nd = 24
  double precision, parameter :: ss(4) = [-1d0,1d0,1d0,-1d0]
  double precision, parameter :: tt(4) = [-1d0,-1d0,1d0,1d0]
  double precision :: gp(3), gw(3), n(4), dn(2,4), b(2,4)
  double precision :: jac(2,2), invj(2,2), detj, geomsize, weight, xi, eta
  double precision :: young, nu, beta, c0, kap, diff, chi, ell, gc, at, kfg
  double precision :: mu, lam, bulk, cg, cold, dg, qg, zg, gradc(2), gradd(2), gradq(2)
  double precision :: ex, ey, es, tr, dev, psi, g, gd, sig0(3), sig(3), peff
  double precision :: aa, ww, hh, fatigue, fatigue_a, bracket, htrial, ah, ih, sh, flux(2)
  double precision :: dex(nd), dey(nd), des(nd), dc(nd), dd(nd), dq(nd), dz(nd)
  double precision :: dgvec(nd), dpsi(nd), da(nd), dh(nd), ds(nd), dpeff(nd)
  double precision :: dsig0(3,nd), dsig(3,nd), df(2,nd), dgradc(2,nd), dgradq(2,nd)
  double precision :: dgradd(2,nd), rr(nd), kk(nd,nd), svnew(27), elastic, zw, source
  integer :: ix, iy, ip, a, j, off, ia, ic, id, iq, iz, lump, mixed_pressure

  rhs = 0d0
  amatrx = 0d0
  energy = 0d0
  lump=1
  mixed_pressure=0
  if (njprop > 0) lump=jprops(1)
  if (njprop > 1) mixed_pressure=jprops(2)
  ! Refuse unsupported sizes, zero time, invalid geometry/data without changing SVARS.
  if (ndofel /= nd .or. nnode /= 4 .or. mcrd < 2 .or. nsvars /= 27 .or. nprops < 11) then
    pnewdt = min(pnewdt,0.5d0)
    return
  end if
  if (.not.ieee_is_finite(dtime) .or. dtime <= 0d0) then
    pnewdt = min(pnewdt,0.5d0)
    return
  end if
  if (.not.all(ieee_is_finite(u)) .or. .not.all(ieee_is_finite(du(:,1))) .or. &
      .not.all(ieee_is_finite(coords)) .or. .not.all(ieee_is_finite(props)) .or. &
      .not.all(ieee_is_finite(svars))) then
    pnewdt = min(pnewdt,0.5d0)
    return
  end if
  young=props(1); nu=props(2); beta=props(3); c0=props(4); kap=props(5)
  diff=props(6); chi=props(7); ell=props(8); gc=props(9); at=props(10); kfg=props(11)
  if (young <= 0d0 .or. nu <= -1d0 .or. nu >= 0.5d0 .or. kap <= 0d0 .or. kap > 1d0 .or. &
      diff < 0d0 .or. ell <= 0d0 .or. gc <= 0d0 .or. at <= 0d0 .or. kfg < 0d0 .or. &
      maxval(abs(u)) > 1d50 .or. any(svars(1:27:3) < 0d0) .or. any(svars(2:27:3) < 0d0)) then
    pnewdt = min(pnewdt,0.5d0)
    return
  end if
  mu = young/(2d0*(1d0+nu))
  lam = young*nu/((1d0+nu)*(1d0-2d0*nu))
  bulk = lam+2d0*mu/3d0
  gp = [-sqrt(0.6d0),0d0,sqrt(0.6d0)]
  gw = [5d0/9d0,8d0/9d0,5d0/9d0]
  geomsize = maxval(coords(1,:))-minval(coords(1,:))
  geomsize = max(geomsize,maxval(coords(2,:))-minval(coords(2,:)))
  rr=0d0; kk=0d0; svnew=svars; elastic=0d0
  ip=0
  do iy=1,3
    eta=gp(iy)
    do ix=1,3
      xi=gp(ix)
      ip=ip+1
      n=0.25d0*(1d0+ss*xi)*(1d0+tt*eta)
      dn(1,:)=0.25d0*ss*(1d0+tt*eta)
      dn(2,:)=0.25d0*tt*(1d0+ss*xi)
      ! jac(i,j)=dx_i/dxi_j, b(i,a)=dN_a/dx_i.
      jac=matmul(coords(1:2,:),transpose(dn))
      detj=jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
      if (.not.ieee_is_finite(detj) .or. detj <= 1d-12*geomsize**2 .or. geomsize <= 0d0) then
        pnewdt=min(pnewdt,0.5d0)
        return
      end if
      invj(1,1)=jac(2,2)/detj; invj(1,2)=-jac(1,2)/detj
      invj(2,1)=-jac(2,1)/detj; invj(2,2)=jac(1,1)/detj
      b=matmul(transpose(invj),dn)
      weight=detj*gw(ix)*gw(iy)
      cg=0d0; cold=0d0; dg=0d0; qg=0d0; zg=0d0
      ex=0d0; ey=0d0; es=0d0; gradc=0d0; gradd=0d0; gradq=0d0
      dex=0d0; dey=0d0; des=0d0; dc=0d0; dd=0d0; dq=0d0; dz=0d0
      dgradc=0d0; dgradq=0d0; dgradd=0d0
      do a=1,4
        off=6*(a-1)
        ex=ex+b(1,a)*u(off+1)
        ey=ey+b(2,a)*u(off+2)
        es=es+0.5d0*(b(2,a)*u(off+1)+b(1,a)*u(off+2))
        cg=cg+n(a)*u(off+3)
        cold=cold+n(a)*(u(off+3)-du(off+3,1))
        dg=dg+n(a)*u(off+4)
        qg=qg+n(a)*u(off+5)
        zg=zg+n(a)*u(off+6)
        gradc=gradc+b(:,a)*u(off+3)
        gradd=gradd+b(:,a)*u(off+4)
        gradq=gradq+b(:,a)*u(off+5)
        dex(off+1)=b(1,a); dey(off+2)=b(2,a)
        des(off+1)=0.5d0*b(2,a); des(off+2)=0.5d0*b(1,a)
        dc(off+3)=n(a); dd(off+4)=n(a); dq(off+5)=n(a); dz(off+6)=n(a)
        dgradc(:,off+3)=b(:,a); dgradq(:,off+5)=b(:,a); dgradd(:,off+4)=b(:,a)
      end do
      ex=ex-beta*(cg-c0); ey=ey-beta*(cg-c0)
      dex=dex-beta*dc; dey=dey-beta*dc
      tr=ex+ey; dev=ex-ey
      sig0=[2d0*mu*ex+lam*tr,2d0*mu*ey+lam*tr,2d0*mu*es]
      dsig0(1,:)=2d0*mu*dex+lam*(dex+dey)
      dsig0(2,:)=2d0*mu*dey+lam*(dex+dey)
      dsig0(3,:)=2d0*mu*des
      peff=(sig0(1)+sig0(2))/3d0
      dpeff=(dsig0(1,:)+dsig0(2,:))/3d0
      g=(1d0-kap)*(1d0-dg)**2+kap
      gd=-2d0*(1d0-kap)*(1d0-dg)
      dgvec=gd*dd
      sig=g*sig0
      do j=1,3
        dsig(j,:)=g*dsig0(j,:)+sig0(j)*dgvec
      end do
      if (mixed_pressure == 1) then
        sig(1)=g*mu*dev+1.5d0*qg
        sig(2)=-g*mu*dev+1.5d0*qg
        dsig(1,:)=mu*(g*(dex-dey)+dev*dgvec)+1.5d0*dq
        dsig(2,:)=-mu*(g*(dex-dey)+dev*dgvec)+1.5d0*dq
      end if
      psi=0.5d0*bulk*max(tr,0d0)**2+mu*(0.5d0*dev**2+2d0*es**2)
      dpsi=bulk*max(tr,0d0)*(dex+dey)+mu*dev*(dex-dey)+4d0*mu*es*des
      ww=g*psi
      off=3*(ip-1)
      aa=svars(off+2)+max(ww-svars(off+3),0d0)
      da=0d0
      if (ww > svars(off+3)) da=g*dpsi+psi*dgvec
      fatigue=1d0; fatigue_a=0d0
      if (aa > at) then
        bracket=max(0d0,1d0-kfg*log10(aa/at))
        fatigue=max(1d-15,bracket**2)
        if (bracket**2 > 1d-15) fatigue_a=-2d0*bracket*kfg/(aa*log(10d0))
      end if
      htrial=2d0*ell*psi/(gc*fatigue)
      hh=max(svars(off+1),htrial)
      dh=0d0
      if (htrial > svars(off+1)) dh=(2d0*ell/(gc*fatigue))*dpsi-htrial*(fatigue_a/fatigue)*da
      svnew(off+1)=hh; svnew(off+2)=aa; svnew(off+3)=ww
      ! Algebraically normalized source constraint. Stable even for H >> 1/eps.
      ! First-order operator split: damage uses accepted H_n.
      ah=(1d0-kap)*svars(off+1)
      dh=0d0
      ih=1d0/(1d0+ah)
      sh=ah*ih
      ds=(1d0-kap)*(ih*ih)*dh
      source=ih*zg+sh*(dg-1d0)
      flux=g*diff*gradc-chi*cg*gradq
      do j=1,nd
        df(:,j)=diff*(g*dgradc(:,j)+dgvec(j)*gradc)-chi*(dc(j)*gradq+cg*dgradq(:,j))
      end do
      do a=1,4
        ia=6*(a-1)+1; ic=ia+2; id=ia+3; iq=ia+4; iz=ia+5
        rr(ia)=rr(ia)+weight*(b(1,a)*sig(1)+b(2,a)*sig(3))
        rr(ia+1)=rr(ia+1)+weight*(b(1,a)*sig(3)+b(2,a)*sig(2))
        rr(ic)=rr(ic)+weight*(n(a)*(cg-cold)/dtime+dot_product(b(:,a),flux))
        rr(iq)=rr(iq)+weight*n(a)*(qg-g*peff)
        kk(ia,:)=kk(ia,:)+weight*(b(1,a)*dsig(1,:)+b(2,a)*dsig(3,:))
        kk(ia+1,:)=kk(ia+1,:)+weight*(b(1,a)*dsig(3,:)+b(2,a)*dsig(2,:))
        kk(ic,:)=kk(ic,:)+weight*(n(a)*dc/dtime+b(1,a)*df(1,:)+b(2,a)*df(2,:))
        kk(iq,:)=kk(iq,:)+weight*n(a)*(dq-g*dpeff-peff*dgvec)
        if (lump == 1) then
          ! Positive nodal quadrature of the zero-order d,z factors; H and
          ! its sensitivities retain their Gauss-point evaluations.
          rr(id)=rr(id)+weight*(ell**2*dot_product(b(:,a),gradd)+n(a)*(u(id)-u(iz)))
          rr(iz)=rr(iz)+weight*n(a)*(ih*u(iz)+sh*(u(id)-1d0))
          kk(id,:)=kk(id,:)+weight*ell**2*(b(1,a)*dgradd(1,:)+b(2,a)*dgradd(2,:))
          kk(id,id)=kk(id,id)+weight*n(a)
          kk(id,iz)=kk(id,iz)-weight*n(a)
          kk(iz,:)=kk(iz,:)+weight*n(a)*(u(id)-1d0-u(iz))*ds
          kk(iz,iz)=kk(iz,iz)+weight*n(a)*ih
          kk(iz,id)=kk(iz,id)+weight*n(a)*sh
        else
          ! Fully integrated alternative, retained for a controlled comparison.
          rr(id)=rr(id)+weight*(ell**2*dot_product(b(:,a),gradd)+n(a)*(dg-zg))
          rr(iz)=rr(iz)+weight*n(a)*source
          kk(id,:)=kk(id,:)+weight*(ell**2*(b(1,a)*dgradd(1,:)+b(2,a)*dgradd(2,:))+n(a)*(dd-dz))
          kk(iz,:)=kk(iz,:)+weight*n(a)*(ih*dz+sh*dd+(dg-1d0-zg)*ds)
        end if
      end do
      zw=mu*(ex**2+ey**2+2d0*es**2)+0.5d0*lam*tr**2
      elastic=elastic+weight*g*zw
    end do
  end do
  if (.not.all(ieee_is_finite(rr)) .or. .not.all(ieee_is_finite(kk)) .or. &
      .not.all(ieee_is_finite(svnew)) .or. .not.ieee_is_finite(elastic)) then
    pnewdt=min(pnewdt,0.5d0)
    return
  end if
  rhs(1:nd,1)=-rr
  amatrx=kk
  svars=svnew
  energy(2)=elastic
end subroutine uel
