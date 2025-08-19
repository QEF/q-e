!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!========================================================================
!                      GGA POTENTIAL DERIVATIVE DRIVERS
!========================================================================
!
!------------------------------------------------------------------------
MODULE qe_drivers_d_gga
  !----------------------------------------------------------------------
  !! Module with QE driver routines that calculates the derivatives of XC
  !! potential.
  !
  USE kind_l,               ONLY: DP
  USE xclib_utils_and_para, ONLY: error_msg, nowarning
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: dgcxc_unpol, dgcxc_spin, d3gcxc
  !
  !
CONTAINS
!
!---------------------------------------------------------------------------
SUBROUTINE dgcxc_unpol( length, r_in, s2_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc )
  !-------------------------------------------------------------------------
  !! This routine computes the derivative of the exchange and correlation
  !! potentials of GGA family.
  !
  USE dft_setting_params,   ONLY: igcx, igcc, is_libxc
  USE qe_drivers_gga,   ONLY: gcxc
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! Number of k-points
  REAL(DP), INTENT(IN), DIMENSION(length) :: r_in
  !! Charge density
  REAL(DP), INTENT(IN), DIMENSION(length) :: s2_in
  !! Square modulus of the density gradient for each k-point
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vrrx
  !! V1x term of the derivative
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vsrx
  !! Cross term for exchange
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vssx
  !! V2x term
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vrrc
  !! V1c term
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vsrc
  !! Cross term for correlation
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vssc
  !! V2c term
  !
  ! ... local variables
  !
  INTEGER :: ir, i1, i2, i3, i4, f1, f2, f3, f4
  INTEGER :: igcx_, igcc_, ierr
  REAL(DP), ALLOCATABLE :: raux(:), s2aux(:), dr(:), s(:), ds(:)
  REAL(DP), ALLOCATABLE :: v1x(:), v2x(:), v1c(:), v2c(:)
  REAL(DP), ALLOCATABLE :: sx(:), sc(:)
  REAL(DP), PARAMETER :: small = 1.E-30_DP
  !
  !$acc data present( r_in, s2_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc )
  !
  ALLOCATE( raux(4*length), s2aux(4*length), dr(length), s(length), ds(length) )
  ALLOCATE( v1x(4*length), v2x(4*length), sx(4*length) )
  ALLOCATE( v1c(4*length), v2c(4*length), sc(4*length) )
  !$acc data create( raux, s2aux, s, dr, ds )
  !$acc data create( v1x, v2x, sx, v1c, v2c, sc )
  !
  igcx_=igcx
  igcc_=igcc
  IF (is_libxc(3)) igcx=0
  IF (is_libxc(4)) igcc=0
  !
  i1 = 1     ;   f1 = length     !4 blocks:  [ rho+dr ,    grho2    ]
  i2 = f1+1  ;   f2 = 2*length   !           [ rho-dr ,    grho2    ]
  i3 = f2+1  ;   f3 = 3*length   !           [ rho    , (grho+ds)^2 ]
  i4 = f3+1  ;   f4 = 4*length   !           [ rho    , (grho-ds)^2 ]
  !
  !$acc parallel loop
  DO ir = 1, length
    s(ir)  = SQRT(s2_in(ir))
    dr(ir) = MIN(1.d-4, 1.d-2*r_in(ir))
    ds(ir) = MIN(1.d-4, 1.d-2*s(ir))
    !
    raux(i1-1+ir) = r_in(ir)+dr(ir)  ;   s2aux(i1-1+ir) = s2_in(ir)
    raux(i2-1+ir) = r_in(ir)-dr(ir)  ;   s2aux(i2-1+ir) = s2_in(ir)
    raux(i3-1+ir) = r_in(ir)         ;   s2aux(i3-1+ir) = (s(ir)+ds(ir))**2
    raux(i4-1+ir) = r_in(ir)         ;   s2aux(i4-1+ir) = (s(ir)-ds(ir))**2
  ENDDO
  !
  CALL gcxc( length*4, raux, s2aux, sx, sc, v1x, v2x, v1c, v2c, ierr )
  !
  !$acc parallel loop
  DO ir = 1, length
    IF ( r_in(ir)>small .AND. s2_in(ir)>small ) THEN
      vrrx(ir) = 0.5_DP * (v1x(i1-1+ir) - v1x(i2-1+ir)) / dr(ir)
      vrrc(ir) = 0.5_DP * (v1c(i1-1+ir) - v1c(i2-1+ir)) / dr(ir)
      !
      vsrx(ir) = 0.25_DP * ((v2x(i1-1+ir) - v2x(i2-1+ir)) / dr(ir) + &
                       (v1x(i3-1+ir) - v1x(i4-1+ir)) / ds(ir) / s(ir))
      vsrc(ir) = 0.25_DP * ((v2c(i1-1+ir) - v2c(i2-1+ir)) / dr(ir) + &
                       (v1c(i3-1+ir) - v1c(i4-1+ir)) / ds(ir) / s(ir))
      !
      vssx(ir) = 0.5_DP * (v2x(i3-1+ir) - v2x(i4-1+ir)) / ds(ir) / s(ir)
      vssc(ir) = 0.5_DP * (v2c(i3-1+ir) - v2c(i4-1+ir)) / ds(ir) / s(ir)
    ELSE
      vrrx(ir) = 0._DP  ;  vrrc(ir) = 0._DP
      vsrx(ir) = 0._DP  ;  vsrc(ir) = 0._DP
      vssx(ir) = 0._DP  ;  vssc(ir) = 0._DP
    ENDIF
  ENDDO
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( raux, s2aux, dr, s, ds )
  DEALLOCATE( v1x, v2x, sx )
  DEALLOCATE( v1c, v2c, sc )
  !
  IF (is_libxc(3)) igcx=igcx_
  IF (is_libxc(4)) igcc=igcc_
  !
  !$acc end data
  !
  IF (ierr/=0 .AND. .NOT.nowarning) CALL xclib_error( 'xc_gcx_', error_msg(ierr), 1 )
  !
  RETURN
  !
END SUBROUTINE dgcxc_unpol
!
!
!--------------------------------------------------------------------------
SUBROUTINE dgcxc_spin( length, r_in, g_in, vrrx, vrsx, vssx, vrrc, vrsc, &
                         vssc, vrzc )
  !------------------------------------------------------------------------
  !! This routine computes the derivative of the exchange and correlation
  !! potentials in the spin-polarized case.
  !
  USE dft_setting_params,   ONLY: igcx, igcc, is_libxc
  USE qe_drivers_gga,   ONLY: gcx_spin, gcc_spin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  !! Number of k-points
  REAL(DP), INTENT(IN), DIMENSION(length,2) :: r_in
  !! Charge density up and down
  REAL(DP), INTENT(IN), DIMENSION(length,3,2) :: g_in
  !! Gradient of the charge density up and down
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vrrx
  !! V1x term of the derivative
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vrsx
  !! Cross term for exchange
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vssx
  !! V2x term of the derivative
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vrrc
  !! V1c term of the derivative
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vrsc
  !! Cross term for correlation
  REAL(DP), INTENT(OUT), DIMENSION(length,2) :: vrzc
  !! Derivative of V1c with respect to zeta
  REAL(DP), INTENT(OUT), DIMENSION(length) :: vssc
  !! V2c term of the derivative
  !
  ! ... local variables
  !
  INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8, ir
  INTEGER :: f1, f2, f3, f4, f5, f6, f7, f8
  ! block delimiters
  INTEGER :: igcx_, igcc_, ierr
  REAL(DP) :: r_up, r_dw, s_up, s_dw, s2_up, s2_dw, rt, zeta, s2t
  REAL(DP) :: dr_up, dr_dw, ds_up, ds_dw, drt, ds, dz
  LOGICAL :: ir_null
  ! used to set output values to zero when input values 
  ! are too small (e.g. rho<eps)
  REAL(DP), ALLOCATABLE :: sx(:), v1x(:,:), v2x(:,:)
  REAL(DP), ALLOCATABLE :: sc(:), v1c(:,:), v2c(:)
  REAL(DP), ALLOCATABLE :: raux(:,:), st(:), s2aux(:,:)
  REAL(DP), ALLOCATABLE :: rtaux(:), s2taux(:), zetaux(:)
  ! auxiliary arrays for gcx- and gcc- routines input
  !
  REAL(DP), PARAMETER :: eps = 1.D-6
  REAL(DP), PARAMETER :: rho_trash = 0.4_DP, zeta_trash = 0.2_DP, &
                         s2_trash = 0.1_DP
  !
  !$acc data present( r_in, g_in, vrrx, vrsx, vssx, vrrc, vrsc, vssc, vrzc )
  !
  igcx_=igcx
  igcc_=igcc
  IF (is_libxc(3)) igcx=0
  IF (is_libxc(4)) igcc=0
  !
  ! ... EXCHANGE
  !
  i1 = 1     ;   f1 = length     !  8 blocks:   [ rho+drup , grho2         ]
  i2 = f1+1  ;   f2 = 2*length   !              [ rho-drup , grho2         ]
  i3 = f2+1  ;   f3 = 3*length   !              [ rho      , (grho+dsup)^2 ]
  i4 = f3+1  ;   f4 = 4*length   !              [ rho      , (grho-dsup)^2 ]
  i5 = f4+1  ;   f5 = 5*length   !              [ rho+drdw , grho2         ]
  i6 = f5+1  ;   f6 = 6*length   !              [ rho-drdw , grho2         ]  
  i7 = f6+1  ;   f7 = 7*length   !              [ rho      , (grho+dsdw)^2 ]
  i8 = f7+1  ;   f8 = 8*length   !              [ rho      , (grho-dsdw)^2 ]    
  !
  ALLOCATE( raux(length*8,2), s2aux(length*8,2) )
  ALLOCATE( sx(length*8), v1x(length*8,2), v2x(length*8,2) )
  !$acc data create( raux, s2aux )
  !$acc data create( sx, v1x, v2x )
  !
  !$acc parallel loop
  DO ir = 1, length
    s2_up = g_in(ir,1,1)**2 + g_in(ir,2,1)**2 + g_in(ir,3,1)**2
    s2_dw = g_in(ir,1,2)**2 + g_in(ir,2,2)**2 + g_in(ir,3,2)**2
    !
    ! ... thresholds
    r_up=r_in(ir,1) ; s_up=SQRT(s2_up)
    IF ( r_in(ir,1)<=eps .OR. SQRT(s2_up)<=eps ) THEN
      r_up=rho_trash ; s2_up=s2_trash ; s_up=SQRT(s2_trash)
    ENDIF
    !
    r_dw=r_in(ir,2) ; s_dw=SQRT(s2_dw)
    IF ( r_in(ir,2)<=eps .OR. SQRT(s2_dw)<=eps ) THEN
      r_dw=rho_trash ; s2_dw=s2_trash ; s_dw=SQRT(s2_trash)
    ENDIF
    !
    dr_up = MIN(1.D-4, 1.D-2*r_up) ; ds_up = MIN(1.D-4, 1.D-2*s_up)
    dr_dw = MIN(1.D-4, 1.D-2*r_dw) ; ds_dw = MIN(1.D-4, 1.D-2*s_dw)
    !
    ! ... up
    raux(i1-1+ir,1) = r_up+dr_up  ;  s2aux(i1-1+ir,1) = s2_up
    raux(i2-1+ir,1) = r_up-dr_up  ;  s2aux(i2-1+ir,1) = s2_up
    raux(i3-1+ir,1) = r_up        ;  s2aux(i3-1+ir,1) = (s_up+ds_up)**2
    raux(i4-1+ir,1) = r_up        ;  s2aux(i4-1+ir,1) = (s_up-ds_up)**2
    !
    raux(i1-1+ir,2) = r_dw        ;  s2aux(i1-1+ir,2) = s2_dw
    raux(i2-1+ir,2) = r_dw        ;  s2aux(i2-1+ir,2) = s2_dw
    raux(i3-1+ir,2) = r_dw        ;  s2aux(i3-1+ir,2) = s2_dw
    raux(i4-1+ir,2) = r_dw        ;  s2aux(i4-1+ir,2) = s2_dw
    ! ... down
    raux(i5-1+ir,1) = r_up        ;  s2aux(i5-1+ir,1) = s2_up
    raux(i6-1+ir,1) = r_up        ;  s2aux(i6-1+ir,1) = s2_up
    raux(i7-1+ir,1) = r_up        ;  s2aux(i7-1+ir,1) = s2_up
    raux(i8-1+ir,1) = r_up        ;  s2aux(i8-1+ir,1) = s2_up
    !
    raux(i5-1+ir,2) = r_dw+dr_dw  ;  s2aux(i5-1+ir,2) = s2_dw
    raux(i6-1+ir,2) = r_dw-dr_dw  ;  s2aux(i6-1+ir,2) = s2_dw
    raux(i7-1+ir,2) = r_dw        ;  s2aux(i7-1+ir,2) = (s_dw+ds_dw)**2
    raux(i8-1+ir,2) = r_dw        ;  s2aux(i8-1+ir,2) = (s_dw-ds_dw)**2
  ENDDO
  !
  CALL gcx_spin( length*8, raux, s2aux, sx, v1x, v2x, ierr )
  !
  !$acc parallel loop
  DO ir = 1, length
    ! ... up
    s2_up = g_in(ir,1,1)**2 + g_in(ir,2,1)**2 + g_in(ir,3,1)**2
    IF ( r_in(ir,1)>eps .AND. SQRT(s2_up)>eps ) THEN
      r_up = r_in(ir,1)
      s_up = SQRT(s2_up)
      dr_up = MIN(1.D-4, 1.D-2*r_up)
      ds_up = MIN(1.D-4, 1.D-2*s_up)
      vrrx(ir,1) = 0.5_DP  *  (v1x(i1-1+ir,1) - v1x(i2-1+ir,1)) / dr_up
      vrsx(ir,1) = 0.25_DP * ((v2x(i1-1+ir,1) - v2x(i2-1+ir,1)) / dr_up + &
                              (v1x(i3-1+ir,1) - v1x(i4-1+ir,1)) / ds_up / s_up)
      vssx(ir,1) = 0.5_DP  *  (v2x(i3-1+ir,1) - v2x(i4-1+ir,1)) / ds_up / s_up
    ELSE
      vrrx(ir,1) = 0._DP
      vrsx(ir,1) = 0._DP
      vssx(ir,1) = 0._DP
    ENDIF
    ! ... down
    s2_dw = g_in(ir,1,2)**2 + g_in(ir,2,2)**2 + g_in(ir,3,2)**2
    IF ( r_in(ir,2)>eps .AND. SQRT(s2_dw)>eps ) THEN
      r_dw = r_in(ir,2)
      s_dw = SQRT(s2_dw)
      dr_dw = MIN(1.D-4, 1.D-2*r_dw)
      ds_dw = MIN(1.D-4, 1.D-2*s_dw)
      vrrx(ir,2) = 0.5_DP  *  (v1x(i5-1+ir,2) - v1x(i6-1+ir,2)) / dr_dw
      vrsx(ir,2) = 0.25_DP * ((v2x(i5-1+ir,2) - v2x(i6-1+ir,2)) / dr_dw + &
                              (v1x(i7-1+ir,2) - v1x(i8-1+ir,2)) / ds_dw / s_dw)
      vssx(ir,2) = 0.5_DP  *  (v2x(i7-1+ir,2) - v2x(i8-1+ir,2)) / ds_dw / s_dw
    ELSE
      vrrx(ir,2) = 0._DP
      vrsx(ir,2) = 0._DP
      vssx(ir,2) = 0._DP
    ENDIF
  ENDDO
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( raux, s2aux  )
  DEALLOCATE( sx, v1x, v2x )
  !
  ! ... CORRELATION
  !
  i1 = 1     ;   f1 = length     !6 blocks: [ rt+dr , s2t       , zeta    ]
  i2 = f1+1  ;   f2 = 2*length   !          [ rt-dr , s2t       , zeta    ]
  i3 = f2+1  ;   f3 = 3*length   !          [ rt    , (st+ds)^2 , zeta    ]
  i4 = f3+1  ;   f4 = 4*length   !          [ rt    , (st-ds)^2 , zeta    ]
  i5 = f4+1  ;   f5 = 5*length   !          [ rt    , grho2     , zeta+dz ]
  i6 = f5+1  ;   f6 = 6*length   !          [ rt    , grho2     , zeta-dz ]  
  !
  ALLOCATE( rtaux(length*6), st(length), s2taux(length*6), zetaux(length*6) )
  ALLOCATE( v1c(length*6,2), v2c(length*6), sc(length*6) )
  !$acc data create( rtaux, st, s2taux, zetaux )
  !$acc data create( v1c, v2c, sc )
  !
  !$acc parallel loop
  DO ir = 1, length
    !
    rt = r_in(ir,1) + r_in(ir,2)
    IF (rt > eps) zeta = (r_in(ir,1)-r_in(ir,2)) / rt
    !
    s2t = (g_in(ir,1,1) + g_in(ir,1,2))**2 + &
          (g_in(ir,2,1) + g_in(ir,2,2))**2 + &
          (g_in(ir,3,1) + g_in(ir,3,2))**2
    st(ir) = SQRT(s2t)
    !
    IF (rt<=eps .OR. ABS(zeta)>1._DP .OR. st(ir)<=eps) THEN
      rt  = rho_trash
      s2t = s2_trash ; st(ir) = SQRT(s2_trash)
      zeta = zeta_trash
    ENDIF
    !
    drt = MIN(1.D-4, 1.D-2 * rt)
    ds = MIN(1.D-4, 1.D-2 * st(ir))
    !dz = MIN(1.D-4, 1.D-2 * ABS(zeta) )
    dz = 1.D-6
    !
    ! ... If zeta is too close to +-1 the derivative is evaluated at a
    ! slightly smaller value.
    zeta = SIGN( MIN(ABS(zeta), (1.0_DP - 2.0_DP*dz)), zeta )
    !
    rtaux(i1-1+ir) = rt+drt ; s2taux(i1-1+ir) = s2t            ; zetaux(i1-1+ir) = zeta
    rtaux(i2-1+ir) = rt-drt ; s2taux(i2-1+ir) = s2t            ; zetaux(i2-1+ir) = zeta
    rtaux(i3-1+ir) = rt     ; s2taux(i3-1+ir) = (st(ir)+ds)**2 ; zetaux(i3-1+ir) = zeta
    rtaux(i4-1+ir) = rt     ; s2taux(i4-1+ir) = (st(ir)-ds)**2 ; zetaux(i4-1+ir) = zeta
    rtaux(i5-1+ir) = rt     ; s2taux(i5-1+ir) = s2t            ; zetaux(i5-1+ir) = zeta+dz
    rtaux(i6-1+ir) = rt     ; s2taux(i6-1+ir) = s2t            ; zetaux(i6-1+ir) = zeta-dz
  ENDDO
  !
  CALL gcc_spin( length*6, rtaux, zetaux, s2taux, sc, v1c, v2c )
  !
  !$acc parallel loop
  DO ir = 1, length
    !
    rt = r_in(ir,1)+r_in(ir,2)
    !
    ir_null = .FALSE.
    zeta = zeta_trash
    IF (rt>=eps) zeta = (r_in(ir,1)-r_in(ir,2)) / rt
    IF (rt<eps .OR. ABS(zeta)>1._DP .OR. st(ir)<eps) ir_null = .TRUE.
    !
    drt = MIN(1.D-4, 1.D-2 * rt)
    ds = MIN(1.D-4, 1.D-2 * st(ir))
    !dz = MIN(1.D-4, 1.D-2 * ABS(zeta) )
    dz = 1.D-6
    !
    IF ( .NOT. ir_null ) THEN
      vrrc(ir,1) = 0.5_DP * (v1c(i1-1+ir,1) - v1c(i2-1+ir,1)) / drt
      vrrc(ir,2) = 0.5_DP * (v1c(i1-1+ir,2) - v1c(i2-1+ir,2)) / drt
      vrsc(ir,1) = 0.5_DP * (v1c(i3-1+ir,1) - v1c(i4-1+ir,1)) / ds/st(ir)
      vrsc(ir,2) = 0.5_DP * (v1c(i3-1+ir,2) - v1c(i4-1+ir,2)) / ds/st(ir)
      vssc(ir)   = 0.5_DP * (v2c(i3-1+ir)   - v2c(i4-1+ir)  ) / ds/st(ir)
      vrzc(ir,1) = 0.5_DP * (v1c(i5-1+ir,1) - v1c(i6-1+ir,1)) / dz
      vrzc(ir,2) = 0.5_DP * (v1c(i5-1+ir,2) - v1c(i6-1+ir,2)) / dz
    ELSE
      vrrc(ir,1) = 0._DP ;  vrrc(ir,2) = 0._DP
      vrsc(ir,1) = 0._DP ;  vrsc(ir,2) = 0._DP
      vrzc(ir,1) = 0._DP ;  vrzc(ir,2) = 0._DP
      vssc(ir)   = 0._DP
    ENDIF
  ENDDO
  !
  !$acc end data
  !$acc end data
  DEALLOCATE( rtaux, st, s2taux, zetaux )
  DEALLOCATE( v1c, v2c, sc )
  !
  IF (is_libxc(3)) igcx=igcx_
  IF (is_libxc(4)) igcc=igcc_
  !
  !$acc end data
  !
  IF (ierr/=0 .AND. .NOT.nowarning) CALL xclib_error( 'xc_gcx_', error_msg(ierr), 2 )
  !
  RETURN
  !
END SUBROUTINE dgcxc_spin
!
!
!-----------------------------------------------------------------------
SUBROUTINE d3gcxc( r, s2, vrrrx, vsrrx, vssrx, vsssx, &
                     vrrrc, vsrrc, vssrc, vsssc )
  !-----------------------------------------------------------------------
  ! wat20101006: 
  !! Calculates all derivatives of the exchange (x) and correlation (c) 
  !! potential in third order of the Exc.
  !
  !    input:       r = rho, s2=|\nabla rho|^2
  !    definition:  E_xc = \int ( f_x(r,s2) + f_c(r,s2) ) dr
  !    output:      vrrrx = d^3(f_x)/d(r)^3
  !                 vsrrx = d^3(f_x)/d(|\nabla r|)d(r)^2 / |\nabla r|
  !                 vssrx = d/d(|\nabla r|) [ &
  !                           d^2(f_x)/d(|\nabla r|)d(r) / |\nabla r| ] &
  !                                                           / |\nabla r|
  !                 vsssx = d/d(|\nabla r|) [ &
  !                           d/d(|\nabla r|) [ &
  !                           d(f_x)/d(|\nabla r|) / |\nabla r| ] &
  !                                                    / |\nabla r| ] &
  !                                                        / |\nabla r|
  !                 same for (c)
  !
  IMPLICIT NONE
  !
  REAL(DP) :: r, s2
  REAL(DP) :: vrrrx, vsrrx, vssrx, vsssx, vrrrc, vsrrc, vssrc, vsssc
  !
  ! ... local variables
  !
  REAL(DP) :: dr, s, ds
  REAL(DP), DIMENSION(4) :: raux, s2aux
  REAL(DP), DIMENSION(4) :: vrrx_rs, vsrx_rs, vssx_rs, vrrc_rs, &
                            vsrc_rs, vssc_rs
  !
  s = SQRT(s2)
  dr = MIN(1.d-4, 1.d-2 * r)
  ds = MIN(1.d-4, 1.d-2 * s)
  !
  raux(1) = r+dr  ; s2aux(1) = s2           !4 blocks:  [ rho+dr ,    grho2    ]
  raux(2) = r-dr  ; s2aux(2) = s2           !           [ rho-dr ,    grho2    ]
  raux(3) = r     ; s2aux(3) = (s+ds)**2    !           [ rho    , (grho+ds)^2 ]
  raux(4) = r     ; s2aux(4) = (s-ds)**2    !           [ rho    , (grho-ds)^2 ]
  !
  CALL dgcxc_unpol( 4, raux, s2aux, vrrx_rs, vsrx_rs, vssx_rs, vrrc_rs, vsrc_rs, vssc_rs )
  !
  vrrrx = 0.5d0  *  (vrrx_rs(1) - vrrx_rs(2)) / dr
  vsrrx = 0.25d0 * ((vsrx_rs(1) - vsrx_rs(2)) / dr &
                   +(vrrx_rs(3) - vrrx_rs(4)) / ds / s)
  vssrx = 0.25d0 * ((vssx_rs(1) - vssx_rs(2)) / dr &
                   +(vsrx_rs(3) - vsrx_rs(4)) / ds / s)
  vsssx = 0.5d0  *  (vssx_rs(3) - vssx_rs(4)) / ds / s
  !
  vrrrc = 0.5d0  *  (vrrc_rs(1) - vrrc_rs(2)) / dr
  vsrrc = 0.25d0 * ((vsrc_rs(1) - vsrc_rs(2)) / dr &
                   +(vrrc_rs(3) - vrrc_rs(4)) / ds / s)
  vssrc = 0.25d0 * ((vssc_rs(1) - vssc_rs(2)) / dr &
                   +(vsrc_rs(3) - vsrc_rs(4)) / ds / s)
  vsssc = 0.5d0  *  (vssc_rs(3) - vssc_rs(4)) / ds / s
  !
  RETURN
  !
END SUBROUTINE d3gcxc
!
END MODULE qe_drivers_d_gga

