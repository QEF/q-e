!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE xc( length, sr_d, sv_d, rho_in, ex_out, ec_out, vx_out, vc_out )
  !-------------------------------------------------------------------------
  !! Wrapper routine. Calls internal XC-driver routines or external ones
  !! from Libxc, depending on the input choice.
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_par_mod,       ONLY: xc_func, xc_info
#endif
  !
  USE kind_l,            ONLY: DP
  USE dft_par_mod,       ONLY: iexch, icorr, is_libxc, rho_threshold_lda, &
                               finite_size_cell_volume_set
  USE qe_drivers_lda_lsda
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: sr_d
  !! spin dimension of rho
  INTEGER, INTENT(IN) :: sv_d
  !! spin dimension of v
  REAL(DP), INTENT(IN) :: rho_in(length,sr_d)
  !! Charge density
  REAL(DP), INTENT(OUT) :: ex_out(length)
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT) :: vx_out(length,sv_d)
  !! \(dE_x(\text{rho})/d\text{rho}  ( NOT d\epsilon_x(\text{rho})/d\text{rho} )
  REAL(DP), INTENT(OUT) :: ec_out(length)
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT) :: vc_out(length,sv_d)
  !! \(dE_c(\text{rho})/d\text{rho}  ( NOT d\epsilon_c(\text{rho})/d\text{rho} )
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  INTEGER :: fkind_x
  REAL(DP) :: amag
  REAL(DP), ALLOCATABLE :: rho_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_lxc(:), vc_lxc(:)
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  !
  REAL(DP), ALLOCATABLE :: arho(:), zeta(:)
  !
  INTEGER  :: ir
  !
  ex_out = 0.0_DP ; vx_out = 0.0_DP
  ec_out = 0.0_DP ; vc_out = 0.0_DP
  !
#if defined(__LIBXC)
  !
  fkind_x = -1
  lengthxc = length
  !
  IF ( ANY(is_libxc(1:2)) ) THEN
    !
    ALLOCATE( rho_lxc(length*sv_d) )
    ALLOCATE( vx_lxc(length*sv_d), vc_lxc(length*sv_d) )
    !
    ! ... set libxc input
    SELECT CASE( sr_d )
    CASE( 1 )
       !
       rho_lxc(:) = ABS(rho_in(:,1))
       !
    CASE( 2 )
       !
       DO ir = 1, length
          rho_lxc(2*ir-1) = (rho_in(ir,1) + rho_in(ir,2)) * 0.5_DP
          rho_lxc(2*ir)   = (rho_in(ir,1) - rho_in(ir,2)) * 0.5_DP
       ENDDO
       !
    CASE( 4 )
       !
       DO ir = 1, length
          amag = SQRT( SUM(rho_in(ir,2:4)**2) )
          rho_lxc(2*ir-1) = (rho_in(ir,1) + amag) * 0.5_DP
          rho_lxc(2*ir)   = (rho_in(ir,1) - amag) * 0.5_DP
       ENDDO
       !
    CASE DEFAULT
       !
       CALL xclib_error( 'xc_LDA', 'Wrong number of spin dimensions', 1 )
       !
    END SELECT
    !
  ENDIF
  !
  !
  ! ... EXCHANGE
  IF ( is_libxc(1) ) THEN
     CALL xc_f03_func_set_dens_threshold( xc_func(1), rho_threshold_lda )
     CALL xc_f03_lda_exc_vxc( xc_func(1), lengthxc, rho_lxc(1), ex_out(1), vx_lxc(1) )
  ENDIF
  !
  ! ... CORRELATION
  IF ( is_libxc(2) ) THEN
     CALL xc_f03_func_set_dens_threshold( xc_func(2), rho_threshold_lda )
     fkind_x = xc_f03_func_info_get_kind( xc_info(2) )
     CALL xc_f03_lda_exc_vxc( xc_func(2), lengthxc, rho_lxc(1), ec_out(1), vc_lxc(1) )
  ENDIF
  !
  IF ( ((.NOT.is_libxc(1)) .OR. (.NOT.is_libxc(2))) &
        .AND. fkind_x/=XC_EXCHANGE_CORRELATION ) THEN
     !
     SELECT CASE( sr_d )
     CASE( 1 )
        !
        IF (iexch==8 .OR. icorr==10) THEN
          IF (.NOT. finite_size_cell_volume_set) CALL xclib_error( 'XC',&
              'finite size corrected exchange used w/o initialization', 1 )
        ENDIF
        CALL xc_lda( length, ABS(rho_in(:,1)), ex_out, ec_out, vx_out(:,1), vc_out(:,1) )
        !
     CASE( 2 )
        !
        ALLOCATE( arho(length), zeta(length) )
        arho = ABS(rho_in(:,1))
        WHERE (arho > rho_threshold_lda) zeta(:) = rho_in(:,2) / arho(:)
        CALL xc_lsda( length, arho, zeta, ex_out, ec_out, vx_out, vc_out )
        DEALLOCATE( arho, zeta )
        !
     CASE( 4 )
        !
        ALLOCATE( arho(length), zeta(length) )
        arho = ABS( rho_in(:,1) )
        WHERE (arho > rho_threshold_lda) zeta(:) = SQRT( rho_in(:,2)**2 + rho_in(:,3)**2 + &
                                             rho_in(:,4)**2 ) / arho(:) ! amag/arho
        CALL xc_lsda( length, arho, zeta, ex_out, ec_out, vx_out, vc_out )
        DEALLOCATE( arho, zeta )
        !
     CASE DEFAULT
        !
        CALL xclib_error( 'xc_LDA', 'Wrong ns input', 2 )
        !
     END SELECT
     !
  ENDIF
  !
  !  ... fill output arrays
  !  
  IF (sv_d == 1) THEN
     IF (is_libxc(1)) vx_out(:,1) = vx_lxc(:)
     IF (is_libxc(2)) vc_out(:,1) = vc_lxc(:)
  ELSE
     IF (is_libxc(1)) THEN
        DO ir = 1, length
           vx_out(ir,1) = vx_lxc(2*ir-1)
           vx_out(ir,2) = vx_lxc(2*ir)
        ENDDO
     ENDIF
     IF (is_libxc(2)) THEN
        DO ir = 1, length
           vc_out(ir,1) = vc_lxc(2*ir-1)
           vc_out(ir,2) = vc_lxc(2*ir)
        ENDDO
     ENDIF
  ENDIF
  !
  IF (ANY(is_libxc(1:2))) THEN
     DEALLOCATE( rho_lxc )
     DEALLOCATE( vx_lxc, vc_lxc )
  ENDIF
  !
#else
  !
  SELECT CASE( sr_d )
  CASE( 1 )
     !
     IF (iexch==8 .OR. icorr==10) THEN
       IF (.NOT. finite_size_cell_volume_set) CALL xclib_error( 'XC',&
           'finite size corrected exchange used w/o initialization', 1 )
     ENDIF
     !
     CALL xc_lda( length, ABS(rho_in(:,1)), ex_out, ec_out, vx_out(:,1), vc_out(:,1) )
     !
  CASE( 2 )
     !
     ALLOCATE( arho(length), zeta(length) )
     !
     arho = ABS(rho_in(:,1))
     WHERE (arho > rho_threshold_lda) zeta(:) = rho_in(:,2) / arho(:)
     !
     CALL xc_lsda( length, arho, zeta, ex_out, ec_out, vx_out, vc_out )
     !
     DEALLOCATE( arho, zeta )
     ! 
  CASE( 4 )
     !
     ALLOCATE( arho(length), zeta(length) )
     !
     arho = ABS( rho_in(:,1) )
     WHERE (arho > rho_threshold_lda) zeta(:) = SQRT( rho_in(:,2)**2 + rho_in(:,3)**2 + &
                                          rho_in(:,4)**2 ) / arho(:) ! amag/arho
     !
     CALL xc_lsda( length, arho, zeta, ex_out, ec_out, vx_out, vc_out )
     !
     DEALLOCATE( arho, zeta )
     !
  CASE DEFAULT
     !
     CALL xclib_error( 'xc_LDA', 'Wrong ns input', 2 )
     !
  END SELECT
  !
#endif
  !
  !
  RETURN
  !
END SUBROUTINE xc
