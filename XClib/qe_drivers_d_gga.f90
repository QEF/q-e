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
  USE kind_l,             ONLY: DP
  USE dft_setting_params, ONLY: igcx, igcc, is_libxc
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
  INTEGER :: i1, i2, i3, i4, f1, f2, f3, f4
  INTEGER :: igcx_, igcc_
  REAL(DP), DIMENSION(length) :: dr, s, ds
  REAL(DP), DIMENSION(4*length) :: raux, s2aux
  REAL(DP), ALLOCATABLE :: v1x(:), v2x(:), v1c(:), v2c(:)
  REAL(DP), ALLOCATABLE :: sx(:), sc(:)
  REAL(DP), PARAMETER :: small = 1.E-30_DP
  !
  ALLOCATE( v1x(4*length), v2x(4*length), sx(4*length) )
  ALLOCATE( v1c(4*length), v2c(4*length), sc(4*length) )
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
  s  = SQRT(s2_in)
  dr = MIN(1.d-4, 1.d-2*r_in)
  ds = MIN(1.d-4, 1.d-2*s)
  !
  raux(i1:f1) = r_in+dr  ;   s2aux(i1:f1) = s2_in
  raux(i2:f2) = r_in-dr  ;   s2aux(i2:f2) = s2_in
  raux(i3:f3) = r_in     ;   s2aux(i3:f3) = (s+ds)**2
  raux(i4:f4) = r_in     ;   s2aux(i4:f4) = (s-ds)**2
  !
  !$acc data copyin( raux, s2aux ) copyout( sx, sc, v1x, v2x, v1c, v2c )
  !$acc host_data use_device( raux, s2aux, sx, sc, v1x, v2x, v1c, v2c )
  CALL gcxc( length*4, raux, s2aux, sx, sc, v1x, v2x, v1c, v2c )
  !$acc end host_data
  !$acc end data
  !
  ! ... to avoid NaN in the next operations
  WHERE( r_in<=small .OR. s2_in<=small )
    dr = 1._DP ; ds = 1._DP ; s = 1._DP
  END WHERE
  !
  vrrx = 0.5_DP * (v1x(i1:f1) - v1x(i2:f2)) / dr
  vrrc = 0.5_DP * (v1c(i1:f1) - v1c(i2:f2)) / dr
  !
  vsrx = 0.25_DP * ((v2x(i1:f1) - v2x(i2:f2)) / dr + &
                    (v1x(i3:f3) - v1x(i4:f4)) / ds / s)
  vsrc = 0.25_DP * ((v2c(i1:f1) - v2c(i2:f2)) / dr + &
                    (v1c(i3:f3) - v1c(i4:f4)) / ds / s)
  !
  vssx = 0.5_DP * (v2x(i3:f3) - v2x(i4:f4)) / ds / s
  vssc = 0.5_DP * (v2c(i3:f3) - v2c(i4:f4)) / ds / s
  !
  DEALLOCATE( v1x, v2x, sx )
  DEALLOCATE( v1c, v2c, sc )
  !
  IF (is_libxc(3)) igcx=igcx_
  IF (is_libxc(4)) igcc=igcc_
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
  INTEGER :: i1, i2, i3, i4, i5, i6, i7, i8
  INTEGER :: f1, f2, f3, f4, f5, f6, f7, f8
  ! block delimiters
  INTEGER :: igcx_, igcc_
  REAL(DP), DIMENSION(length,2) :: r, s, s2
  REAL(DP), DIMENSION(length,2) :: drup, drdw, dsup, dsdw
  REAL(DP), ALLOCATABLE :: sx(:), v1x(:,:), v2x(:,:)
  REAL(DP), ALLOCATABLE :: sc(:), v1c(:,:), v2c(:)
  REAL(DP), DIMENSION(length) :: rt, zeta, st, s2t
  REAL(DP), DIMENSION(length) :: dr, ds, dz
  REAL(DP), DIMENSION(length,2) :: null_v
  ! used to set output values to zero when input values 
  ! are too small (e.g. rho<eps)
  ! 
  REAL(DP), ALLOCATABLE :: raux(:,:), s2aux(:,:)
  REAL(DP), ALLOCATABLE :: rtaux(:), s2taux(:), zetaux(:)
  ! auxiliary arrays for gcx- and gcc- routines input
  !
  REAL(DP), PARAMETER :: eps = 1.D-6
  REAL(DP), PARAMETER :: rho_trash = 0.4_DP, zeta_trash = 0.2_DP, &
                         s2_trash = 0.1_DP
  !
  igcx_=igcx
  igcc_=igcc
  IF (is_libxc(3)) igcx=0
  IF (is_libxc(4)) igcc=0
  !
  vrrx = 0.0_DP ; vrsx = 0.0_DP ; vssx = 0.0_DP
  vrrc = 0.0_DP ; vrsc = 0.0_DP ; vrzc = 0.0_DP
  vssc = 0.0_DP
  !
  ! ... EXCHANGE
  !
  i1 = 1     ;   f1 = length     !8 blocks(x2): [ rho+drup , grho2         ]
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
  !
  s2(:,1) = g_in(:,1,1)**2 + g_in(:,2,1)**2 + g_in(:,3,1)**2 !up
  s2(:,2) = g_in(:,1,2)**2 + g_in(:,2,2)**2 + g_in(:,3,2)**2 !down
  s = SQRT(s2)
  r = r_in
  !
  null_v = 1.0_DP
  !
  ! ... thresholds
  WHERE ( r_in(:,1)<=eps .OR. s(:,1)<=eps )
     r(:,1) = rho_trash
     s2(:,1) = s2_trash ; s(:,1) = SQRT(s2_trash)
     null_v(:,1) = 0.0_DP
  END WHERE
  !
  WHERE ( r_in(:,2)<=eps .OR. s(:,2)<=eps )
     r(:,2) = rho_trash
     s2(:,2) = s2_trash ; s(:,2) = SQRT(s2_trash)
     null_v(:,2) = 0.0_DP
  END WHERE
  !
  drup = 0.0_DP ;  drdw = 0.0_DP
  dsup = 0.0_DP ;  dsdw = 0.0_DP
  !
  drup(:,1) = MIN(1.D-4, 1.D-2*r(:,1)) ; dsup(:,1) = MIN(1.D-4, 1.D-2*s(:,1))
  drdw(:,2) = MIN(1.D-4, 1.D-2*r(:,2)) ; dsdw(:,2) = MIN(1.D-4, 1.D-2*s(:,2))
  !
  ! ... up
  raux(i1:f1,:) = r+drup ;  s2aux(i1:f1,:) = s2
  raux(i2:f2,:) = r-drup ;  s2aux(i2:f2,:) = s2
  raux(i3:f3,:) = r      ;  s2aux(i3:f3,:) = (s+dsup)**2
  raux(i4:f4,:) = r      ;  s2aux(i4:f4,:) = (s-dsup)**2
  ! ... down
  raux(i5:f5,:) = r+drdw ;  s2aux(i5:f5,:) = s2
  raux(i6:f6,:) = r-drdw ;  s2aux(i6:f6,:) = s2
  raux(i7:f7,:) = r      ;  s2aux(i7:f7,:) = (s+dsdw)**2
  raux(i8:f8,:) = r      ;  s2aux(i8:f8,:) = (s-dsdw)**2
  !
  !
  !$acc data copyin( raux, s2aux ) copyout( sx, v1x, v2x )
  !$acc host_data use_device( raux, s2aux, sx, v1x, v2x )
  CALL gcx_spin( length*8, raux, s2aux, sx, v1x, v2x )
  !$acc end host_data
  !$acc end data
  !
  ! ... up
  vrrx(:,1) = 0.5_DP  *  (v1x(i1:f1,1) - v1x(i2:f2,1)) / drup(:,1)
  vrsx(:,1) = 0.25_DP * ((v2x(i1:f1,1) - v2x(i2:f2,1)) / drup(:,1) + &
                         (v1x(i3:f3,1) - v1x(i4:f4,1)) / dsup(:,1) / s(:,1))
  vssx(:,1) = 0.5_DP  *  (v2x(i3:f3,1) - v2x(i4:f4,1)) / dsup(:,1) / s(:,1)
  ! ... down
  vrrx(:,2) = 0.5_DP  *  (v1x(i5:f5,2) - v1x(i6:f6,2)) / drdw(:,2)
  vrsx(:,2) = 0.25_DP * ((v2x(i5:f5,2) - v2x(i6:f6,2)) / drdw(:,2) + &
                                     (v1x(i7:f7,2) - v1x(i8:f8,2)) / dsdw(:,2) / s(:,2))
  vssx(:,2) = 0.5_DP  *  (v2x(i7:f7,2) - v2x(i8:f8,2)) / dsdw(:,2) / s(:,2)
  !
  vrrx(:,1) = vrrx(:,1)*null_v(:,1) ;  vrrx(:,2) = vrrx(:,2)*null_v(:,2)
  vrsx(:,1) = vrsx(:,1)*null_v(:,1) ;  vrsx(:,2) = vrsx(:,2)*null_v(:,2)
  vssx(:,1) = vssx(:,1)*null_v(:,1) ;  vssx(:,2) = vssx(:,2)*null_v(:,2)
  !
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
  ALLOCATE( rtaux(length*6), s2taux(length*6), zetaux(length*6) )
  ALLOCATE( v1c(length*6,2), v2c(length*6), sc(length*6) )
  !
  rt(:) = r_in(:,1) + r_in(:,2)
  !
  null_v = 1.0_DP
  !
  WHERE (rt > eps)
     zeta = (r_in(:,1) - r_in(:,2)) / rt(:)
  ELSEWHERE
     zeta = zeta_trash
     null_v(:,1) = 0.0_DP
  END WHERE
  !
  s2t = (g_in(:,1,1) + g_in(:,1,2))**2 + &
        (g_in(:,2,1) + g_in(:,2,2))**2 + &
        (g_in(:,3,1) + g_in(:,3,2))**2
  st = SQRT(s2t)
  !
  WHERE (rt<eps .OR. ABS(zeta)>1._DP .OR. st<eps)
     rt(:)  = rho_trash
     s2t(:) = s2_trash ; st = SQRT(s2_trash)
     zeta(:) = zeta_trash
     null_v(:,1) = 0.0_DP
  END WHERE
  !
  dr = MIN(1.D-4, 1.D-2 * rt)
  ds = MIN(1.D-4, 1.D-2 * st)
  !dz = MIN(1.D-4, 1.D-2 * ABS(zeta) )
  dz = 1.D-6
  !
  ! ... If zeta is too close to +-1 the derivative is evaluated at a
  ! slightly smaller value.
  zeta = SIGN( MIN(ABS(zeta), (1.0_DP - 2.0_DP*dz)), zeta )
  !
  rtaux(i1:f1) = rt+dr ;  s2taux(i1:f1) = s2t        ;  zetaux(i1:f1) = zeta
  rtaux(i2:f2) = rt-dr ;  s2taux(i2:f2) = s2t        ;  zetaux(i2:f2) = zeta
  rtaux(i3:f3) = rt    ;  s2taux(i3:f3) = (st+ds)**2 ;  zetaux(i3:f3) = zeta
  rtaux(i4:f4) = rt    ;  s2taux(i4:f4) = (st-ds)**2 ;  zetaux(i4:f4) = zeta
  rtaux(i5:f5) = rt    ;  s2taux(i5:f5) = s2t        ;  zetaux(i5:f5) = zeta+dz
  rtaux(i6:f6) = rt    ;  s2taux(i6:f6) = s2t        ;  zetaux(i6:f6) = zeta-dz
  !
  !$acc data copyin( rtaux, zetaux, s2taux ) copyout( sc, v1c, v2c )
  !$acc host_data use_device( rtaux, zetaux, s2taux, sc, v1c, v2c )
  CALL gcc_spin( length*6, rtaux, zetaux, s2taux, sc, v1c, v2c )
  !$acc end host_data
  !$acc end data
  !
  vrrc(:,1) = 0.5_DP * (v1c(i1:f1,1) - v1c(i2:f2,1)) / dr    * null_v(:,1)
  vrrc(:,2) = 0.5_DP * (v1c(i1:f1,2) - v1c(i2:f2,2)) / dr    * null_v(:,1)
  vrsc(:,1) = 0.5_DP * (v1c(i3:f3,1) - v1c(i4:f4,1)) / ds/st * null_v(:,1)
  vrsc(:,2) = 0.5_DP * (v1c(i3:f3,2) - v1c(i4:f4,2)) / ds/st * null_v(:,1)
  vssc(:)   = 0.5_DP * (v2c(i3:f3)   - v2c(i4:f4)  ) / ds/st * null_v(:,1)
  vrzc(:,1) = 0.5_DP * (v1c(i5:f5,1) - v1c(i6:f6,1)) / dz    * null_v(:,1)
  vrzc(:,2) = 0.5_DP * (v1c(i5:f5,2) - v1c(i6:f6,2)) / dz    * null_v(:,1)
  !
  IF (is_libxc(3)) igcx=igcx_
  IF (is_libxc(4)) igcc=igcc_
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

