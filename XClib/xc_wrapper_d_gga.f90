
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE dgcxc( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss, gpu_args_ )
  !---------------------------------------------------------------------
  !! Wrapper routine. Calls dgcx-driver routines from internal libraries
  !! or from the external libxc, depending on the input choice.
  !
  USE kind_l, ONLY: DP
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
  LOGICAL, OPTIONAL, INTENT(IN) :: gpu_args_
  !! whether you wish to run on gpu in case use_gpu is true
  !
  LOGICAL :: gpu_args
  !
  gpu_args = .FALSE.
  IF ( PRESENT(gpu_args_) ) gpu_args = gpu_args_
  !
  IF ( gpu_args ) THEN
    !
    !$acc data present( r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
    CALL dgcxc_( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
    !$acc end data
    !
  ELSE
    !
    !$acc data copyin( r_in, g_in ), copyout( dvxc_rr, dvxc_sr, dvxc_ss )
    CALL dgcxc_( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
    !$acc end data
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE
!
!---------------------------------------------------------------------
SUBROUTINE dgcxc_( length, sp, r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
  !---------------------------------------------------------------------
  !! Wrapper routine. Calls dgcx-driver routines from internal libraries
  !! or from the external libxc, depending on the input choice.
  !
  USE constants_l,          ONLY: e2
  USE kind_l,               ONLY: DP
  USE dft_setting_params,   ONLY: igcx, igcc, is_libxc, rho_threshold_gga, &
                                  grho_threshold_gga, rho_threshold_lda,   &
                                  ishybrid, exx_started, exx_fraction
  USE qe_drivers_d_gga
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_setting_params,   ONLY: xc_func, xc_info
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
  INTEGER :: fkind
#if defined(__LIBXC)
  REAL(DP), ALLOCATABLE :: rho_lxc(:)
  REAL(DP), ALLOCATABLE :: v2rho2_x(:), v2rhosigma_x(:), v2sigma2_x(:)
  REAL(DP), ALLOCATABLE :: v2rho2_c(:), v2rhosigma_c(:), v2sigma2_c(:)
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  !
  LOGICAL :: fkind_is_XC
  INTEGER :: k, length_dlxc
  REAL(DP) :: rht, zeta, xcoef
  REAL(DP), ALLOCATABLE :: sigma(:)
  REAL(DP), PARAMETER :: small=1.E-10_DP, rho_trash=0.5_DP
  REAL(DP), PARAMETER :: epsr=1.0d-6, epsg=1.0d-6
  !
  !$acc data present( r_in, g_in, dvxc_rr, dvxc_sr, dvxc_ss )
  !
  IF ( ANY(.NOT.is_libxc(3:4)) ) THEN
    rho_threshold_gga = small ;  grho_threshold_gga = small
  ENDIF
  !
  !$acc kernels
  dvxc_rr(:,:,:) = 0._DP
  dvxc_sr(:,:,:) = 0._DP
  dvxc_ss(:,:,:) = 0._DP
  !$acc end kernels
  !
  fkind = -1
  fkind_is_XC = .FALSE.
  !
#if defined(__LIBXC)
  !
  IF ( ANY(is_libxc(3:4)) ) THEN
    !
    lengthxc = length
    !
    length_dlxc = length
    IF (sp == 2) length_dlxc = length*3
    !
    ALLOCATE( rho_lxc(length*sp), sigma(length_dlxc) )
    !$acc enter data create( rho_lxc, sigma )
    !
    ! ... set libxc input
    SELECT CASE( sp )
    CASE( 1 )
      !
      !$acc parallel loop
      DO k = 1, length
        rho_lxc(k) = r_in(k,1)
        sigma(k) = g_in(k,1,1)**2 + g_in(k,2,1)**2 + g_in(k,3,1)**2
      ENDDO
      !
    CASE( 2 )
      !
      !$acc parallel loop
      DO k = 1, length
        rho_lxc(2*k-1) = r_in(k,1)
        rho_lxc(2*k)   = r_in(k,2)
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
    !$acc update self(rho_lxc, sigma)
    !
  ENDIF
  !
  ! ... LIBXC DERIVATIVE FOR EXCHANGE
  !
  IF ( is_libxc(3) .AND. igcx/=0 ) THEN
    ALLOCATE( v2rho2_x(length_dlxc), v2rhosigma_x(length_dlxc*sp), &
              v2sigma2_x(length_dlxc*sp) )
    v2rho2_x = 0._DP ;  v2rhosigma_x = 0._DP ;  v2sigma2_x = 0._DP
    CALL xc_f03_func_set_dens_threshold( xc_func(3), epsr )
    CALL xc_f03_gga_fxc( xc_func(3), lengthxc, rho_lxc(1), sigma(1), v2rho2_x(1), &
                         v2rhosigma_x(1), v2sigma2_x(1) )
    !$acc data copyin( v2rho2_x, v2rhosigma_x, v2sigma2_x )
    !
    xcoef = 1.d0
    IF ( ishybrid .AND. exx_started .AND. exx_fraction>0.d0) xcoef = 1.d0-exx_fraction
    !
    IF (sp==1) THEN
      !$acc parallel loop
      DO k = 1, length
        IF ( rho_lxc(k)>small .AND. SQRT(ABS(sigma(k)))>small ) THEN
          IF ( rho_lxc(k)>rho_threshold_lda ) THEN
            dvxc_rr(k,1,1) = xcoef * e2 * v2rho2_x(k)
            dvxc_sr(k,1,1) = xcoef * e2 * v2rhosigma_x(k)*2._DP
          ENDIF
          dvxc_ss(k,1,1) = xcoef * e2 * v2sigma2_x(k)*4._DP
        ENDIF
      ENDDO
    ELSEIF (sp==2) THEN
      !$acc parallel loop
      DO k = 1, length
        IF ( (r_in(k,1)>epsr .AND. SQRT(ABS(sigma(3*k-2)))>epsg) .AND. &
             (r_in(k,2)>epsr .AND. SQRT(ABS(sigma(3*k)))  >epsg) ) THEN
          dvxc_rr(k,1,1) = xcoef * e2 * v2rho2_x(3*k-2)
          dvxc_ss(k,1,1) = xcoef * e2 * v2sigma2_x(6*k-5)*4._DP
          dvxc_rr(k,2,2) = xcoef * e2 * v2rho2_x(3*k)
          dvxc_ss(k,2,2) = xcoef * e2 * v2sigma2_x(6*k)*4._DP
          dvxc_rr(k,1,2) = xcoef * e2 * v2rho2_x(3*k-1)
          dvxc_sr(k,1,1) = xcoef * e2 * v2rhosigma_x(6*k-5)*2._DP
          dvxc_rr(k,2,1) = xcoef * e2 * v2rho2_x(3*k-1)
          dvxc_sr(k,2,2) = xcoef * e2 * v2rhosigma_x(6*k)*2._DP
        ENDIF
      ENDDO
    ENDIF
    !$acc end data
    DEALLOCATE( v2rho2_x, v2rhosigma_x, v2sigma2_x )
    !
  ENDIF
  !
  ! ... LIBXC DERIVATIVE FOR CORRELATION
  !
  IF ( is_libxc(4) .AND. igcc/=0 ) THEN
    ALLOCATE( v2rho2_c(length_dlxc), v2rhosigma_c(length_dlxc*sp), &
              v2sigma2_c(length_dlxc*sp) )
    ! ... DERIVATIVE FOR CORRELATION
    v2rho2_c = 0._DP ;  v2rhosigma_c = 0._DP ;  v2sigma2_c = 0._DP
    fkind = xc_f03_func_info_get_kind( xc_info(4) )
    CALL xc_f03_func_set_dens_threshold( xc_func(4), epsr )
    CALL xc_f03_gga_fxc( xc_func(4), lengthxc, rho_lxc(1), sigma(1), v2rho2_c(1), &
                           v2rhosigma_c(1), v2sigma2_c(1) )
    !$acc data copyin( v2rho2_c, v2rhosigma_c, v2sigma2_c )
    !
    IF (sp==1) THEN
      !$acc parallel loop
      DO k = 1, length
        IF ( rho_lxc(k)>small .AND. SQRT(ABS(sigma(k)))>small ) THEN
          IF ( rho_lxc(k)>rho_threshold_lda ) THEN
            dvxc_rr(k,1,1) = dvxc_rr(k,1,1) + e2 * v2rho2_c(k)
            dvxc_sr(k,1,1) = dvxc_sr(k,1,1) + e2 * v2rhosigma_c(k)*2._DP
          ENDIF
          dvxc_ss(k,1,1) = dvxc_ss(k,1,1) + e2 * v2sigma2_c(k)*4._DP
        ENDIF
      ENDDO
    ELSEIF (sp==2) THEN
      !$acc parallel loop
      DO k = 1, length
        IF ( (r_in(k,1)>epsr .AND. SQRT(ABS(sigma(3*k-2)))>epsg) .AND. &
             (r_in(k,2)>epsr .AND. SQRT(ABS(sigma(3*k)))  >epsg) ) THEN
          dvxc_rr(k,1,1) = dvxc_rr(k,1,1) + e2 * v2rho2_c(3*k-2)
          dvxc_rr(k,1,2) = dvxc_rr(k,1,2) + e2 * v2rho2_c(3*k-1)
          dvxc_rr(k,2,1) = dvxc_rr(k,2,1) + e2 * v2rho2_c(3*k-1)
          dvxc_rr(k,2,2) = dvxc_rr(k,2,2) + e2 * v2rho2_c(3*k)
          dvxc_sr(k,1,1) = dvxc_sr(k,1,1) + e2 * v2rhosigma_c(6*k-5)*2.d0
          dvxc_ss(k,1,1) = dvxc_ss(k,1,1) + e2 * v2sigma2_c(6*k)*4.d0
          dvxc_sr(k,1,2) = dvxc_sr(k,1,2) + e2 * v2rhosigma_c(6*k-4)
          dvxc_sr(k,2,1) = dvxc_sr(k,2,1) + e2 * v2rhosigma_c(6*k-1)
          dvxc_ss(k,1,2) = dvxc_ss(k,1,2) + e2 * v2sigma2_c(6*k-2)
          dvxc_ss(k,2,1) = dvxc_ss(k,2,1) + e2 * v2sigma2_c(6*k-2)
          dvxc_sr(k,2,2) = dvxc_sr(k,2,2) + e2 * v2rhosigma_c(6*k)*2.d0
          dvxc_ss(k,2,2) = dvxc_ss(k,2,2) + e2 * v2sigma2_c(6*k)*4.d0
        ENDIF
      ENDDO
    ENDIF
    !
    !$acc end data
    DEALLOCATE( v2rho2_c, v2rhosigma_c, v2sigma2_c )
    !
  ENDIF
  !
  fkind_is_XC = (fkind==XC_EXCHANGE_CORRELATION)
  !
#endif
  !
  ! ... QE DERIVATIVE FOR EXCHANGE AND CORRELATION
  !
  IF ( ((.NOT.is_libxc(3).AND.igcx/=0) .OR. (.NOT.is_libxc(4).AND.igcc/=0)) &
        .AND. (.NOT.fkind_is_XC) ) THEN
    !
    ALLOCATE( vrrx(length,sp), vsrx(length,sp), vssx(length,sp) )
    ALLOCATE( vrrc(length,sp), vsrc(length,sp), vssc(length) )
    !$acc data create( vrrx, vsrx, vssx, vrrc, vsrc, vssc )
    !
    IF ( sp == 1 ) THEN
       !
       IF (.NOT. ALLOCATED(sigma)) THEN
         ALLOCATE( sigma(length) )
         !$acc enter data create(sigma)
       ENDIF
       !
       !$acc parallel loop
       DO k = 1, length
         sigma(k) = g_in(k,1,1)**2 + g_in(k,2,1)**2 + g_in(k,3,1)**2
       ENDDO
       !
       CALL dgcxc_unpol( length, r_in(:,1), sigma, vrrx(:,1), vsrx(:,1), vssx(:,1), &
                         vrrc(:,1), vsrc(:,1), vssc )
       !
       !$acc parallel loop
       DO k = 1, length
         dvxc_rr(k,1,1) = dvxc_rr(k,1,1) + e2 * (vrrx(k,1) + vrrc(k,1))
         dvxc_sr(k,1,1) = dvxc_sr(k,1,1) + e2 * (vsrx(k,1) + vsrc(k,1))
         dvxc_ss(k,1,1) = dvxc_ss(k,1,1) + e2 * (vssx(k,1) + vssc(k)  )
       ENDDO
       !
    ELSEIF ( sp == 2 ) THEN
       !
       ALLOCATE( vrzc(length,sp) )
       !$acc data create( vrzc )
       !
       CALL dgcxc_spin( length, r_in, g_in, vrrx, vsrx, vssx, vrrc, vsrc, vssc, vrzc )
       !
       !$acc parallel loop
       DO k = 1, length
         rht = r_in(k,1) + r_in(k,2)
         IF (rht > epsr) THEN
           zeta = (r_in(k,1) - r_in(k,2))/rht
           !
           dvxc_rr(k,1,1) = dvxc_rr(k,1,1) + e2*(vrrx(k,1) + vrrc(k,1) + &
                                                 vrzc(k,1)*(1.d0 - zeta)/rht)
           dvxc_rr(k,1,2) = dvxc_rr(k,1,2) + e2*(vrrc(k,1) - vrzc(k,1)*(1.d0 + zeta)/rht)
           dvxc_rr(k,2,1) = dvxc_rr(k,2,1) + e2*(vrrc(k,2) + vrzc(k,2)*(1.d0 - zeta)/rht)
           dvxc_rr(k,2,2) = dvxc_rr(k,2,2) + e2*(vrrx(k,2) + vrrc(k,2) - &
                                                 vrzc(k,2)*(1.d0 + zeta)/rht)
         ENDIF
         !
         dvxc_sr(k,1,1) = dvxc_sr(k,1,1) + e2 * (vsrx(k,1) + vsrc(k,1))
         dvxc_sr(k,1,2) = dvxc_sr(k,1,2) + e2 * vsrc(k,1)
         dvxc_sr(k,2,1) = dvxc_sr(k,2,1) + e2 * vsrc(k,2)
         dvxc_sr(k,2,2) = dvxc_sr(k,2,2) + e2 * (vsrx(k,2) + vsrc(k,2))
         !
         dvxc_ss(k,1,1) = dvxc_ss(k,1,1) + e2 * (vssx(k,1) + vssc(k))
         dvxc_ss(k,1,2) = dvxc_ss(k,1,2) + e2 * vssc(k)
         dvxc_ss(k,2,1) = dvxc_ss(k,2,1) + e2 * vssc(k)
         dvxc_ss(k,2,2) = dvxc_ss(k,2,2) + e2 * (vssx(k,2) + vssc(k))
       ENDDO
       !
       !$acc end data
       DEALLOCATE( vrzc )
       !
    ENDIF
    !
    !$acc end data
    DEALLOCATE( vrrx, vsrx, vssx )
    DEALLOCATE( vrrc, vsrc, vssc )
    !
  ENDIF
  !
#if defined(__LIBXC)
  IF ( ANY(is_libxc(3:4))) THEN
    !$acc exit data delete(rho_lxc)
    DEALLOCATE( rho_lxc )
  ENDIF
#endif
  IF ( ALLOCATED(sigma) ) THEN
    !$acc exit data delete(sigma)
    DEALLOCATE( sigma )
  ENDIF
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE dgcxc_
