!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
!------- DRIVERS FOR DERIVATIVES OF XC POTENTIAL (GGA CASE) ------------
!-----------------------------------------------------------------------
!
!---------------------------------------------------------------------
SUBROUTINE dgcxc( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
  !---------------------------------------------------------------------
  !! Wrapper routine. Calls dgcx-driver routines from internal libraries
  !! or from the external libxc, depending on the input choice.
  !
  USE constants_l,      ONLY: e2
  USE kind_l,           ONLY: DP
  USE dft_par_mod,      ONLY: igcx, igcc, is_libxc, rho_threshold_gga, &
                              grho_threshold_gga
  USE qe_drivers_d_gga
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_par_mod,      ONLY: xc_func, xc_info
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: sp
  !! number of spin components
  REAL(DP), INTENT(IN) :: r_in(length,sp)
  !! charge density
  REAL(DP), INTENT(IN) :: g_in(length,3,sp)
  !! gradient
  REAL(DP), INTENT(OUT) :: dvxc_rr(length,sp,sp), dvxc_sr(length,sp,sp), &
                           dvxc_ss(length,sp,sp)
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: vrrx(:,:), vsrx(:,:), vssx(:,:)
  REAL(DP), ALLOCATABLE :: vrrc(:,:), vsrc(:,:), vssc(:), vrzc(:,:) 
  !
#if defined(__LIBXC)
  INTEGER :: fkind
  REAL(DP), ALLOCATABLE :: rho_lbxc(:)
  REAL(DP), ALLOCATABLE :: v2rho2_x(:), v2rhosigma_x(:), v2sigma2_x(:) 
  REAL(DP), ALLOCATABLE :: v2rho2_c(:), v2rhosigma_c(:), v2sigma2_c(:) 
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  !
  INTEGER :: k, ir, length_lxc, length_dlxc
  REAL(DP) :: rht, zeta
  REAL(DP), ALLOCATABLE :: sigma(:)
  REAL(DP), PARAMETER :: small = 1.E-10_DP, rho_trash = 0.5_DP
  REAL(DP), PARAMETER :: epsr=1.0d-6, epsg=1.0d-10
  !
  IF ( ANY(.NOT.is_libxc(3:4)) ) THEN
    rho_threshold_gga = small ;  grho_threshold_gga = small
  ENDIF
  !
#if defined(__LIBXC)
  !
  IF ( ANY(is_libxc(3:4)) ) THEN
    !
    lengthxc = length
    !
    length_lxc = length*sp
    length_dlxc = length
    IF (sp == 2) length_dlxc = length*3
    !
    ALLOCATE( rho_lbxc(length_lxc), sigma(length_dlxc) )
    !
    ! ... set libxc input
    SELECT CASE( sp )
    CASE( 1 )
      !
      DO k = 1, length
        rho_lbxc(k) = r_in(k,1) 
        sigma(k) = g_in(k,1,1)**2 + g_in(k,2,1)**2 + g_in(k,3,1)**2
      ENDDO
      !
    CASE( 2 )
      !
      DO k = 1, length
        rho_lbxc(2*k-1) = r_in(k,1)
        rho_lbxc(2*k)   = r_in(k,2)
        !
        sigma(3*k-2) = g_in(k,1,1)**2 + g_in(k,2,1)**2 + g_in(k,3,1)**2
        sigma(3*k-1) = g_in(k,1,1) * g_in(k,1,2) + g_in(k,2,1) * g_in(k,2,2) + &
                       g_in(k,3,1) * g_in(k,3,2)
        sigma(3*k)   = g_in(k,1,2)**2 + g_in(k,2,2)**2 + g_in(k,3,2)**2
      ENDDO
      !
    CASE( 4 )
      !
      CALL xclib_error( 'dgcxc', 'The Libxc derivative of the XC potential &
                            &is not available for noncollinear case', 1 )
      !
    CASE DEFAULT
      !
      CALL xclib_error( 'dgcxc', 'Wrong number of spin dimensions', 2 )
      !
    END SELECT
    !
  ENDIF
  !
  IF ( is_libxc(3) ) THEN
    ALLOCATE( v2rho2_x(length_dlxc), v2rhosigma_x(length_dlxc*sp), &
              v2sigma2_x(length_dlxc*sp) )
    ! ... DERIVATIVE FOR EXCHANGE
    v2rho2_x = 0._DP ;  v2rhosigma_x = 0._DP ;  v2sigma2_x = 0._DP
    IF (igcx /= 0) THEN 
      CALL xc_f03_gga_fxc( xc_func(3), lengthxc, rho_lbxc(1), sigma(1), v2rho2_x(1), &
                           v2rhosigma_x(1), v2sigma2_x(1) )
    ENDIF
  ENDIF
  !
  IF ( is_libxc(4) ) THEN
    ALLOCATE( v2rho2_c(length_dlxc), v2rhosigma_c(length_dlxc*sp), &
              v2sigma2_c(length_dlxc*sp) )
    ! ... DERIVATIVE FOR CORRELATION
    v2rho2_c = 0._DP ;  v2rhosigma_c = 0._DP ;  v2sigma2_c = 0._DP
    IF (igcc /= 0) THEN 
      fkind = xc_f03_func_info_get_kind( xc_info(4) )
      CALL xc_f03_gga_fxc( xc_func(4), lengthxc, rho_lbxc(1), sigma(1), v2rho2_c(1), &
                           v2rhosigma_c(1), v2sigma2_c(1) )
    ENDIF
  ENDIF
  !
  IF (ANY(is_libxc(3:4)))  DEALLOCATE( rho_lbxc, sigma )
  !
  dvxc_rr = 0._DP
  dvxc_sr = 0._DP
  dvxc_ss = 0._DP
  !
  IF ( ((.NOT.is_libxc(3)) .OR. (.NOT.is_libxc(4))) &
        .AND. fkind/=XC_EXCHANGE_CORRELATION ) THEN
    !
    ALLOCATE( vrrx(length,sp), vsrx(length,sp), vssx(length,sp) )
    ALLOCATE( vrrc(length,sp), vsrc(length,sp), vssc(length) )
    !
    IF ( sp == 1 ) THEN
       !
       ALLOCATE( sigma(length) )
       sigma(:) = g_in(:,1,1)**2 + g_in(:,2,1)**2 + g_in(:,3,1)**2
       !
       CALL dgcxc_unpol( length, r_in(:,1), sigma, vrrx(:,1), vsrx(:,1), vssx(:,1), &
                         vrrc(:,1), vsrc(:,1), vssc )
       !
       DEALLOCATE( sigma )
       !
    ELSEIF ( sp == 2 ) THEN
       !
       ALLOCATE( vrzc(length,sp) )
       !
       CALL dgcxc_spin( length, r_in, g_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc, vrzc )
       !
    ENDIF
    !
  ENDIF
  !
  IF ( sp == 1 ) THEN
    !
    IF ( ((.NOT.is_libxc(3)) .OR. (.NOT.is_libxc(4))) ) THEN
      dvxc_rr(:,1,1) = e2 * (vrrx(:,1) + vrrc(:,1))
      dvxc_sr(:,1,1) = e2 * (vsrx(:,1) + vsrc(:,1))
      dvxc_ss(:,1,1) = e2 * (vssx(:,1) + vssc(:)  )
      !
      DEALLOCATE( vrrx, vsrx, vssx )
      DEALLOCATE( vrrc, vsrc, vssc )
      !
    ENDIF  
    !
    IF ( is_libxc(3) ) THEN
      dvxc_rr(:,1,1) = dvxc_rr(:,1,1) + e2 * v2rho2_x(:)
      dvxc_sr(:,1,1) = dvxc_sr(:,1,1) + e2 * v2rhosigma_x(:)*2._DP
      dvxc_ss(:,1,1) = dvxc_ss(:,1,1) + e2 * v2sigma2_x(:)*4._DP
    ENDIF
    !
    IF ( is_libxc(4) ) THEN
      dvxc_rr(:,1,1) = dvxc_rr(:,1,1) + e2 * v2rho2_c(:)
      dvxc_sr(:,1,1) = dvxc_sr(:,1,1) + e2 * v2rhosigma_c(:)*2._DP
      dvxc_ss(:,1,1) = dvxc_ss(:,1,1) + e2 * v2sigma2_c(:)*4._DP
    ENDIF
    !
  ELSEIF ( sp == 2 ) THEN
    !
    IF ( ((.NOT.is_libxc(3)) .OR. (.NOT.is_libxc(4))) ) THEN
      !
      DO k = 1, length
        rht = r_in(k,1) + r_in(k,2)
        IF (rht > epsr) THEN
          zeta = (r_in(k,1) - r_in(k,2)) / rht
          !
          dvxc_rr(k,1,1) = e2*(vrrx(k,1) + vrrc(k,1) + vrzc(k,1)*(1._DP - zeta)/rht)
          dvxc_rr(k,1,2) = e2*(vrrc(k,1) - vrzc(k,1)*(1._DP + zeta)/rht)
          dvxc_rr(k,2,1) = e2*(vrrc(k,2) + vrzc(k,2)*(1._DP - zeta)/rht)
          dvxc_rr(k,2,2) = e2*(vrrx(k,2) + vrrc(k,2) - vrzc(k,2)*(1._DP + zeta)/rht)
        ENDIF
        !
        dvxc_sr(k,1,1) = e2*(vsrx(k,1) + vsrc(k,1))
        dvxc_sr(k,1,2) = e2*vsrc(k,1)
        dvxc_sr(k,2,1) = e2*vsrc(k,2)
        dvxc_sr(k,2,2) = e2*(vsrx(k,2) + vsrc(k,2))
        !
        dvxc_ss(k,1,1) = e2*(vssx(k,1) + vssc(k))
        dvxc_ss(k,1,2) = e2*vssc(k)
        dvxc_ss(k,2,1) = e2*vssc(k)
        dvxc_ss(k,2,2) = e2*(vssx(k,2) + vssc(k))
        !
      ENDDO 
      !
      DEALLOCATE( vrrx, vsrx, vssx )
      DEALLOCATE( vrrc, vsrc, vssc )
      DEALLOCATE( vrzc )
      !
    ENDIF
    !   
    IF ( is_libxc(3) ) THEN   
      !   
      DO k = 1, length   
        rht = r_in(k,1) + r_in(k,2)   
        IF (rht > epsr) THEN   
          dvxc_rr(k,1,1) = dvxc_rr(k,1,1) + e2 * v2rho2_x(3*k-2)   
          dvxc_rr(k,1,2) = dvxc_rr(k,1,2) + e2 * v2rho2_x(3*k-1)   
          dvxc_rr(k,2,1) = dvxc_rr(k,2,1) + e2 * v2rho2_x(3*k-1)   
          dvxc_rr(k,2,2) = dvxc_rr(k,2,2) + e2 * v2rho2_x(3*k)   
        ENDIF   
        dvxc_sr(k,1,1) = dvxc_sr(k,1,1) + e2 * v2rhosigma_x(6*k-5)*2._DP
        dvxc_ss(k,1,1) = dvxc_ss(k,1,1) + e2 * v2sigma2_x(6*k-5)*4._DP
        dvxc_sr(k,2,2) = dvxc_sr(k,2,2) + e2 * v2rhosigma_x(6*k)*2._DP  
        dvxc_ss(k,2,2) = dvxc_ss(k,2,2) + e2 * v2sigma2_x(6*k)*4._DP
      ENDDO   
      !   
      DEALLOCATE( v2rho2_x, v2rhosigma_x, v2sigma2_x )   
      !   
    ENDIF   
    !   
    !   
    IF ( is_libxc(4) ) THEN   
      !   
      DO k = 1, length   
        rht = r_in(k,1) + r_in(k,2)   
        IF (rht > epsr) THEN   
          dvxc_rr(k,1,1) = dvxc_rr(k,1,1) + e2 * v2rho2_c(3*k-2)   
          dvxc_rr(k,1,2) = dvxc_rr(k,1,2) + e2 * v2rho2_c(3*k-1)   
          dvxc_rr(k,2,1) = dvxc_rr(k,2,1) + e2 * v2rho2_c(3*k-1)   
          dvxc_rr(k,2,2) = dvxc_rr(k,2,2) + e2 * v2rho2_c(3*k)   
        ENDIF   
        dvxc_sr(k,1,1) = dvxc_sr(k,1,1) + e2 * v2rhosigma_c(6*k-5)*2.d0   
        dvxc_ss(k,1,1) = dvxc_ss(k,1,1) + e2 * v2sigma2_c(6*k)*4.d0   
        dvxc_sr(k,1,2) = dvxc_sr(k,1,2) + e2 * v2rhosigma_c(6*k-4)   
        dvxc_sr(k,2,1) = dvxc_sr(k,2,1) + e2 * v2rhosigma_c(6*k-1)   
        dvxc_ss(k,1,2) = dvxc_ss(k,1,2) + e2 * v2sigma2_c(6*k-2)   
        dvxc_ss(k,2,1) = dvxc_ss(k,2,1) + e2 * v2sigma2_c(6*k-2)   
        dvxc_sr(k,2,2) = dvxc_sr(k,2,1) + e2 * v2rhosigma_c(6*k)*2.d0   
        dvxc_ss(k,2,2) = dvxc_ss(k,2,2) + e2 * v2sigma2_c(6*k)*4.d0   
      ENDDO   
      !   
      DEALLOCATE( v2rho2_c, v2rhosigma_c, v2sigma2_c )   
      !   
    ENDIF   
    !   
  ENDIF
  !
#else
  !
  ALLOCATE( vrrx(length,sp), vsrx(length,sp), vssx(length,sp) )
  ALLOCATE( vrrc(length,sp), vsrc(length,sp), vssc(length) )
  !
  SELECT CASE( sp )
  CASE( 1 )
     !
     ALLOCATE( sigma(length) )
     sigma(:) = g_in(:,1,1)**2 + g_in(:,2,1)**2 + g_in(:,3,1)**2
     CALL dgcxc_unpol( length, r_in(:,1), sigma, vrrx(:,1), vsrx(:,1), vssx(:,1), &
                       vrrc(:,1), vsrc(:,1), vssc )
     DEALLOCATE( sigma )
     !
     dvxc_rr(:,1,1) = e2*(vrrx(:,1) + vrrc(:,1))
     dvxc_sr(:,1,1) = e2*(vsrx(:,1) + vsrc(:,1))
     dvxc_ss(:,1,1) = e2*(vssx(:,1) + vssc(:)  )     
     !
  CASE( 2 )
     !
     ALLOCATE( vrzc(length,sp) )
     !
     CALL dgcxc_spin( length, r_in, g_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc, vrzc )
     !
     DO k = 1, length
        !
        rht = r_in(k,1) + r_in(k,2)
        IF (rht > epsr) THEN
           zeta = (r_in(k,1) - r_in(k,2))/rht
           !
           dvxc_rr(k,1,1) = e2*(vrrx(k,1) + vrrc(k,1) + vrzc(k,1)*(1.d0 - zeta)/rht)
           dvxc_rr(k,1,2) = e2*(vrrc(k,1) - vrzc(k,1)*(1.d0 + zeta)/rht)
           dvxc_rr(k,2,1) = e2*(vrrc(k,2) + vrzc(k,2)*(1.d0 - zeta)/rht)
           dvxc_rr(k,2,2) = e2*(vrrx(k,2) + vrrc(k,2) - vrzc(k,2)*(1.d0 + zeta)/rht)
        ENDIF
        !
        dvxc_sr(k,1,1) = e2 * (vsrx(k,1) + vsrc(k,1))   
        dvxc_sr(k,1,2) = e2 * vsrc(k,1)   
        dvxc_sr(k,2,1) = e2 * vsrc(k,2)   
        dvxc_sr(k,2,2) = e2 * (vsrx(k,2) + vsrc(k,2))   
        !
        dvxc_ss(k,1,1) = e2 * (vssx(k,1) + vssc(k))   
        dvxc_ss(k,1,2) = e2 * vssc(k)   
        dvxc_ss(k,2,1) = e2 * vssc(k)   
        dvxc_ss(k,2,2) = e2 * (vssx(k,2) + vssc(k))
     ENDDO
     !
     DEALLOCATE( vrzc )
     !
  CASE DEFAULT
     !
     CALL xclib_error( 'dgcxc', 'Wrong ns input', 4 )
     !
  END SELECT
  !
  DEALLOCATE( vrrx, vsrx, vssx )
  DEALLOCATE( vrrc, vsrc, vssc )
  !
#endif
  !
  !
  RETURN
  !
END SUBROUTINE dgcxc
