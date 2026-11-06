! Independent higher-order comparison: affine Q2 x Q2, four physical fields.
! Direct broken Hessian of Q2 displacement in grad(g*p); no q,z auxiliaries.
! Reference only. The delivered/scored element remains four-node Q1.
subroutine uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords,mcrd,nnode, &
 u,du,v,acc,jtype,time,dtime,kstep,kinc,jelem,params,ndload,jdltyp,adlmag,predef,npredf, &
 lflags,mlvarx,ddlmag,mdload,pnewdt,jprops,njprop,period)
 use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
 implicit none
 integer,intent(in)::ndofel,nrhs,nsvars,nprops,mcrd,nnode,jtype,kstep,kinc,jelem,ndload,npredf,mlvarx,mdload,njprop
 integer,intent(in)::jdltyp(mdload,*),lflags(*),jprops(*)
 double precision,intent(in)::props(nprops),coords(mcrd,nnode),u(ndofel),du(mlvarx,*)
 double precision,intent(in)::v(ndofel),acc(ndofel),time(2),dtime,params(*),adlmag(mdload,*),predef(2,npredf,nnode),ddlmag(mdload,*),period
 double precision,intent(inout)::svars(nsvars),pnewdt
 double precision,intent(out)::rhs(mlvarx,nrhs),amatrx(ndofel,ndofel),energy(8)
 integer,parameter::nd=36
 double precision,parameter::gp(5)=[-.906179845938664d0,-.538469310105683d0,0d0,.538469310105683d0,.906179845938664d0]
 double precision,parameter::gw(5)=[.236926885056189d0,.478628670499366d0,.568888888888889d0,.478628670499366d0,.236926885056189d0]
 double precision::lx(3),ly(3),dx(3),dy(3),ddp(3),n(9),b(2,9),hs(3,9)
 double precision::E,nu,beta,c0,kap,D,chi,ell,gc,at,kfg,mu,lam,bulk,pc,hx,hy,xi,eta,w
 double precision::ex,ey,es,tr,dev,c,cp,dmg,g,gd,p,pg(2),dg(2),cg(2),qg(2),flux(2),sig0(3),sig(3)
 double precision::psi,Wcur,alpha,f,ht,aa,elastic
 double precision::vx(nd),vy(nd),vs(nd),vc(nd),vd(nd),dgvec(nd),vp(nd),vpg(2,nd),vcg(2,nd),vdg(2,nd)
 double precision::vqg(2,nd),vf(2,nd),vsg(3,nd),r(nd),k(nd,nd),sn(75)
 integer::ix,iy,a,ax,ay,off,ip,j,iu,ic,id
 rhs=0d0; amatrx=0d0; energy=0d0
 if(ndofel/=36.or.nnode/=9.or.nsvars/=75.or.nprops<11.or.dtime<=0d0)then
   pnewdt=.5d0; return
 end if
 E=props(1);nu=props(2);beta=props(3);c0=props(4);kap=props(5);D=props(6)
 chi=props(7);ell=props(8);gc=props(9);at=props(10);kfg=props(11)
 mu=E/(2*(1+nu));lam=E*nu/((1+nu)*(1-2*nu));bulk=lam+2*mu/3
 pc=2*(mu+lam)/3
 hx=maxval(coords(1,:))-minval(coords(1,:));hy=maxval(coords(2,:))-minval(coords(2,:))
 if(hx<=0d0.or.hy<=0d0)then
   pnewdt=.5d0;return
 end if
 ddp=[1d0,-2d0,1d0];r=0d0;k=0d0;sn=svars;elastic=0d0;ip=0
 do iy=1,5
  eta=gp(iy)
  ly=[eta*(eta-1)/2,1-eta**2,eta*(eta+1)/2];dy=[eta-.5d0,-2*eta,eta+.5d0]
  do ix=1,5
   xi=gp(ix);ip=ip+1
   lx=[xi*(xi-1)/2,1-xi**2,xi*(xi+1)/2];dx=[xi-.5d0,-2*xi,xi+.5d0]
   w=hx*hy*gw(ix)*gw(iy)/4
   a=0
   do ay=1,3
    do ax=1,3
     a=a+1;n(a)=lx(ax)*ly(ay)
     b(:,a)=[dx(ax)*ly(ay)*2/hx,lx(ax)*dy(ay)*2/hy]
     hs(:,a)=[ddp(ax)*ly(ay)*4/hx**2,dx(ax)*dy(ay)*4/(hx*hy),lx(ax)*ddp(ay)*4/hy**2]
    end do
   end do
   ex=0d0;ey=0d0;es=0d0;c=0d0;cp=0d0;dmg=0d0;cg=0d0;dg=0d0;pg=0d0
   vx=0d0;vy=0d0;vs=0d0;vc=0d0;vd=0d0;vpg=0d0;vcg=0d0;vdg=0d0
   do a=1,9
    off=4*(a-1)
    ex=ex+b(1,a)*u(off+1);ey=ey+b(2,a)*u(off+2)
    es=es+(b(2,a)*u(off+1)+b(1,a)*u(off+2))/2
    c=c+n(a)*u(off+3);cp=cp+n(a)*(u(off+3)-du(off+3,1));dmg=dmg+n(a)*u(off+4)
    cg=cg+b(:,a)*u(off+3);dg=dg+b(:,a)*u(off+4)
    pg=pg+[hs(1,a)*u(off+1)+hs(2,a)*u(off+2),hs(2,a)*u(off+1)+hs(3,a)*u(off+2)]
    vx(off+1)=b(1,a);vy(off+2)=b(2,a);vs(off+1)=b(2,a)/2;vs(off+2)=b(1,a)/2
    vc(off+3)=n(a);vd(off+4)=n(a);vcg(:,off+3)=b(:,a);vdg(:,off+4)=b(:,a)
    vpg(:,off+1)=pc*hs(1:2,a);vpg(:,off+2)=pc*hs(2:3,a);vpg(:,off+3)=-2*pc*beta*b(:,a)
   end do
   ex=ex-beta*(c-c0);ey=ey-beta*(c-c0);vx=vx-beta*vc;vy=vy-beta*vc
   tr=ex+ey;dev=ex-ey;p=pc*tr;vp=pc*(vx+vy);pg=pc*(pg-2*beta*cg)
   g=(1-kap)*(1-dmg)**2+kap;gd=-2*(1-kap)*(1-dmg);dgvec=gd*vd
   sig0=[2*mu*ex+lam*tr,2*mu*ey+lam*tr,2*mu*es];sig=g*sig0
   vsg(1,:)=g*(2*mu*vx+lam*(vx+vy))+sig0(1)*dgvec
   vsg(2,:)=g*(2*mu*vy+lam*(vx+vy))+sig0(2)*dgvec
   vsg(3,:)=g*2*mu*vs+sig0(3)*dgvec
   qg=g*pg+p*gd*dg
   do j=1,nd
    vqg(:,j)=dgvec(j)*pg+g*vpg(:,j)+vp(j)*gd*dg+p*2*(1-kap)*vd(j)*dg+p*gd*vdg(:,j)
   end do
   flux=g*D*cg-chi*c*qg
   do j=1,nd
    vf(:,j)=D*(dgvec(j)*cg+g*vcg(:,j))-chi*(vc(j)*qg+c*vqg(:,j))
   end do
   psi=bulk*max(tr,0d0)**2/2+mu*(dev**2/2+2*es**2)
   Wcur=g*psi;off=3*(ip-1);alpha=svars(off+2)+max(Wcur-svars(off+3),0d0)
   f=1d0
   if(alpha>at)f=max(1d-15,max(0d0,1-kfg*log10(alpha/at))**2)
   ht=max(svars(off+1),2*ell*psi/(gc*f))
   sn(off+1)=ht;sn(off+2)=alpha;sn(off+3)=Wcur
   ! The same first-order history split as the candidate, but a different
   ! spatial method, direct Q2 Hessians and the original consistent d equation.
   aa=(1-kap)*svars(off+1)
   do a=1,9
    iu=4*(a-1)+1;ic=iu+2;id=iu+3
    r(iu)=r(iu)+w*(b(1,a)*sig(1)+b(2,a)*sig(3))
    r(iu+1)=r(iu+1)+w*(b(1,a)*sig(3)+b(2,a)*sig(2))
    r(ic)=r(ic)+w*(n(a)*(c-cp)/dtime+dot_product(b(:,a),flux))
    r(id)=r(id)+w*(ell**2*dot_product(b(:,a),dg)+n(a)*(dmg-aa*(1-dmg)))
    k(iu,:)=k(iu,:)+w*(b(1,a)*vsg(1,:)+b(2,a)*vsg(3,:))
    k(iu+1,:)=k(iu+1,:)+w*(b(1,a)*vsg(3,:)+b(2,a)*vsg(2,:))
    k(ic,:)=k(ic,:)+w*(n(a)*vc/dtime+b(1,a)*vf(1,:)+b(2,a)*vf(2,:))
    k(id,:)=k(id,:)+w*(ell**2*(b(1,a)*vdg(1,:)+b(2,a)*vdg(2,:))+n(a)*(1+aa)*vd)
   end do
   elastic=elastic+w*g*(mu*(ex**2+ey**2+2*es**2)+lam*tr**2/2)
  end do
 end do
 if(.not.all(ieee_is_finite(r)).or..not.all(ieee_is_finite(k)).or..not.all(ieee_is_finite(sn)))then
  pnewdt=.5d0;return
 end if
 rhs(1:nd,1)=-r;amatrx=k;svars=sn;energy(2)=elastic
end subroutine uel
