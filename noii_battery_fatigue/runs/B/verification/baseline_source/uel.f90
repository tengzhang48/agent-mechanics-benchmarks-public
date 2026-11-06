! Checked standard Abaqus/Standard UEL entry. All physics/assembly is generated
! in generated_model.for. This wrapper owns no persistent memory.
subroutine UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, PROPS, NPROPS, &
               COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME, DTIME, KSTEP, KINC, &
               JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, PREDEF, NPREDF, LFLAGS, &
               MLVARX, DDLMAG, MDLOAD, PNEWDT, JPROPS, NJPROP, PERIOD)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE, KSTEP, KINC
  integer :: JELEM, NDLOAD, NPREDF, MLVARX, MDLOAD, NJPROP
  integer :: JDLTYP(MDLOAD,*), LFLAGS(*), JPROPS(*)
  double precision :: RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL), SVARS(NSVARS), ENERGY(8)
  double precision :: PROPS(NPROPS), COORDS(MCRD,NNODE), U(NDOFEL), DU(MLVARX,*)
  double precision :: V(NDOFEL), A(NDOFEL), TIME(2), DTIME, PARAMS(3)
  double precision :: ADLMAG(MDLOAD,*), PREDEF(2,NPREDF,NNODE), DDLMAG(MDLOAD,*)
  double precision :: PNEWDT, PERIOD
  double precision :: trial_state(12), deriv(4,2), shp(4), jac(2,2), detj, xi, eta
  double precision :: geom_scale, p_trial
  integer :: i, j, k, l

  RHS(1:MLVARX,1:NRHS) = 0.0d0
  AMATRX = 0.0d0
  ENERGY = 0.0d0
  ! These procedures are outside this quasistatic, single-RHS contract.
  if (NRHS /= 1) return
  if (LFLAGS(1) /= 1 .and. LFLAGS(1) /= 2 .and. LFLAGS(1) /= 71 .and. &
      LFLAGS(1) /= 72 .and. LFLAGS(1) /= 73) return
  if (LFLAGS(4) /= 0) return
  if (LFLAGS(3) /= 1 .and. LFLAGS(3) /= 2 .and. LFLAGS(3) /= 5) return

  ! Fail closed without writing SVARS. DTIME=0 is a diagnostic, not an update.
  if (NDOFEL /= 20 .or. NNODE /= 4 .or. MCRD < 2 .or. NSVARS /= 12 .or. &
      NPROPS /= 11 .or. MLVARX < NDOFEL) goto 900
  if (.not. ieee_is_finite(DTIME) .or. DTIME <= 0.0d0) goto 900
  if (.not. all(ieee_is_finite(PROPS))) goto 900
  if (PROPS(1) <= 0.0d0 .or. PROPS(2) <= -1.0d0 .or. PROPS(2) >= 0.5d0) goto 900
  if (PROPS(5) <= 0.0d0 .or. PROPS(5) > 1.0d0 .or. PROPS(6) < 0.0d0) goto 900
  if (PROPS(8) <= 0.0d0 .or. PROPS(9) <= 0.0d0 .or. PROPS(10) <= 0.0d0) goto 900
  if (.not. all(ieee_is_finite(U))) goto 900
  if (.not. all(ieee_is_finite(DU(1:NDOFEL,1)))) goto 900
  if (.not. all(ieee_is_finite(COORDS))) goto 900
  if (.not. all(ieee_is_finite(SVARS))) goto 900
  if (.not. all(ieee_is_finite(TIME))) goto 900
  if (any(SVARS < 0.0d0)) goto 900
  if (maxval(abs(U)) > 1.0d50 .or. maxval(abs(PROPS)) > 1.0d50) goto 900

  ! A bilinear Q4 Jacobian determinant is affine in the natural coordinates;
  ! positive determinants at every corner imply positivity on the whole quad.
  geom_scale = maxval((COORDS(1:2,2)-COORDS(1:2,1))**2) &
             + maxval((COORDS(1:2,4)-COORDS(1:2,1))**2)
  do j = 1, 2
    eta = dble(2*j-3)
    do i = 1, 2
      xi = dble(2*i-3)
      call shape_quad4(xi, eta, shp, deriv)
      jac = 0.0d0
      do k = 1, 4
        do l = 1, 2
          jac(l,:) = jac(l,:) + deriv(k,l)*COORDS(1:2,k)
        end do
      end do
      detj = jac(1,1)*jac(2,2)-jac(1,2)*jac(2,1)
      if (.not. ieee_is_finite(detj) .or. detj <= 1.0d-12*geom_scale) goto 900
    end do
  end do

  trial_state = SVARS
  p_trial = PNEWDT
  call CHEM_UEL_CORE(RHS, AMATRX, trial_state, ENERGY, NDOFEL, NRHS, NSVARS, &
                    PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, &
                    TIME, DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, &
                    ADLMAG, PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, &
                    p_trial, JPROPS, NJPROP, PERIOD)
  if (.not. ieee_is_finite(p_trial) .or. p_trial < 1.0d0) goto 900
  if (.not. all(ieee_is_finite(RHS(1:NDOFEL,1)))) goto 900
  if (.not. all(ieee_is_finite(AMATRX))) goto 900
  if (.not. all(ieee_is_finite(trial_state))) goto 900
  do i = 1, 12, 3
    if (trial_state(i) < SVARS(i) .or. trial_state(i+1) < SVARS(i+1)) goto 900
  end do
  if (LFLAGS(3) == 1) SVARS = trial_state
  return

900 continue
  RHS(1:MLVARX,1:NRHS) = 0.0d0
  AMATRX = 0.0d0
  ENERGY = 0.0d0
  PNEWDT = 0.5d0
end subroutine UEL
