! Handwritten UEL guard. Generated bulk/history and accepted_damage together
! supply the residual and consistent tangent for the first-order time split.
subroutine UEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, &
  PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME, &
  DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, &
  PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, &
  JPROPS, NJPROP, PERIOD)
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: NDOFEL, NRHS, NSVARS, NPROPS, MCRD, NNODE, JTYPE
  integer :: KSTEP, KINC, JELEM, NDLOAD, NPREDF, MLVARX, MDLOAD, NJPROP
  integer :: JDLTYP(MDLOAD,*), LFLAGS(*), JPROPS(*)
  double precision :: RHS(MLVARX,*), AMATRX(NDOFEL,NDOFEL), SVARS(NSVARS), ENERGY(8)
  double precision :: PROPS(*), COORDS(MCRD,NNODE), U(NDOFEL), DU(MLVARX,*)
  double precision :: V(NDOFEL), A(NDOFEL), TIME(2), DTIME, PARAMS(*)
  double precision :: ADLMAG(MDLOAD,*), PREDEF(2,NPREDF,NNODE), DDLMAG(MDLOAD,*)
  double precision :: PNEWDT, PERIOD
  double precision :: saved(NSVARS)
  integer :: i
  logical :: ok

  RHS(1:MLVARX,1:NRHS) = 0d0
  AMATRX = 0d0
  ENERGY = 0d0
  ! Zero-time calls have no meaningful concentration rate or fatigue update.
  if (.not. ieee_is_finite(DTIME)) then
    PNEWDT = 0.5d0
    return
  endif
  if (DTIME <= 0d0) then
    PNEWDT = 0.5d0
    return
  endif
  if (NDOFEL /= 24 .or. NNODE /= 4 .or. MCRD /= 2 .or. &
      NPROPS /= 11 .or. NSVARS < 12 .or. NRHS /= 1 .or. MLVARX < NDOFEL) then
    PNEWDT = 0.5d0
    return
  endif
  if ((LFLAGS(1) /= 1 .and. LFLAGS(1) /= 2 .and. LFLAGS(1) /= 71 .and. &
       LFLAGS(1) /= 72 .and. LFLAGS(1) /= 73) .or. LFLAGS(4) /= 0 .or. &
      (LFLAGS(3) /= 1 .and. LFLAGS(3) /= 2 .and. LFLAGS(3) /= 5)) then
    PNEWDT = 0.5d0
    return
  endif
  ok = all(ieee_is_finite(PROPS(1:NPROPS))) .and. &
       all(ieee_is_finite(COORDS)) .and. all(ieee_is_finite(U)) .and. &
       all(ieee_is_finite(DU(1:NDOFEL,1))) .and. all(ieee_is_finite(SVARS))
  if (.not. ok) then
    PNEWDT = 0.5d0
    return
  endif
  ok = PROPS(1) > 0d0 .and. PROPS(2) > -1d0 .and. PROPS(2) < 0.5d0 .and. &
       PROPS(5) > 0d0 .and. PROPS(5) <= 1d0 .and. PROPS(6) >= 0d0 .and. &
       PROPS(8) > 0d0 .and. PROPS(9) > 0d0 .and. PROPS(10) > 0d0 .and. PROPS(11) >= 0d0
  if (.not. ok) then
    PNEWDT = 0.5d0
    return
  endif
  saved = SVARS
  if (any(SVARS(1:12:3) < 0d0) .or. any(SVARS(2:12:3) < 0d0)) then
    PNEWDT = 0.5d0
    return
  endif
  call CHEMO_KERNEL(RHS, AMATRX, SVARS, ENERGY, NDOFEL, NRHS, NSVARS, &
    PROPS, NPROPS, COORDS, MCRD, NNODE, U, DU, V, A, JTYPE, TIME, &
    DTIME, KSTEP, KINC, JELEM, PARAMS, NDLOAD, JDLTYP, ADLMAG, &
    PREDEF, NPREDF, LFLAGS, MLVARX, DDLMAG, MDLOAD, PNEWDT, &
    JPROPS, NJPROP, PERIOD)
  call accepted_damage(RHS,AMATRX,COORDS,U,saved,PROPS,NDOFEL,MLVARX,PNEWDT)
  ok = all(ieee_is_finite(RHS(1:NDOFEL,1))) .and. &
       all(ieee_is_finite(AMATRX)) .and. all(ieee_is_finite(SVARS))
  if (.not. ok .or. PNEWDT < 1d0) then
    SVARS = saved
    RHS(1:MLVARX,1:NRHS) = 0d0
    AMATRX = 0d0
    ENERGY = 0d0
    PNEWDT = 0.5d0
    return
  endif
  ! Guard the declared irreversible quantities against invalid old states.
  do i = 1, 4
    if (SVARS(3*i-2) < saved(3*i-2) .or. SVARS(3*i-1) < saved(3*i-1)) then
      SVARS = saved
      RHS(1:MLVARX,1:NRHS) = 0d0
      AMATRX = 0d0
      ENERGY = 0d0
      PNEWDT = 0.5d0
      return
    endif
  enddo
end subroutine UEL
