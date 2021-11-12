!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------------------
SUBROUTINE xc( length, srd, svd, rho_in, ex_out, ec_out, vx_out, vc_out, run_on_gpu_ )
  !--------------------------------------------------------------------------------------
  !! Wrapper routine to \(\texttt{xc_}\) or \(\texttt{xc_gpu}\).
  !
  USE kind_l, ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: srd
  !! spin dimension of rho
  INTEGER, INTENT(IN) :: svd
  !! spin dimension of v
  REAL(DP), INTENT(IN) :: rho_in(length,srd)
  !! Charge density
  REAL(DP), INTENT(OUT) :: ex_out(length)
  !! \(\epsilon_x(rho)\) ( NOT \(E_x(\text{rho})\) )
  REAL(DP), INTENT(OUT) :: vx_out(length,svd)
  !! \(dE_x(\text{rho})/d\text{rho}  ( NOT d\epsilon_x(\text{rho})/d\text{rho} )
  REAL(DP), INTENT(OUT) :: ec_out(length)
  !! \(\epsilon_c(rho)\) ( NOT \(E_c(\text{rho})\) )
  REAL(DP), INTENT(OUT) :: vc_out(length,svd)
  !! \(dE_c(\text{rho})/d\text{rho}  ( NOT d\epsilon_c(\text{rho})/d\text{rho} )
  LOGICAL, OPTIONAL, INTENT(IN) :: run_on_gpu_
  !! whether you wish to run on gpu in case use_gpu is true
  !
  LOGICAL :: run_on_gpu
  !
  run_on_gpu = .FALSE.
  IF ( PRESENT(run_on_gpu_) ) run_on_gpu = run_on_gpu_
  !
  IF ( run_on_gpu ) THEN
    !
    !$acc data deviceptr( rho_in(length,srd), ex_out(length), ec_out(length), &
    !$acc&                vx_out(length,svd), vc_out(length,svd) )
    CALL xc_gpu( length, srd, svd, rho_in, ex_out, ec_out, vx_out, vc_out )
    !$acc end data
    !
  ELSE
    !
    CALL xc_( length, srd, svd, rho_in, ex_out, ec_out, vx_out, vc_out )
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE
!
!
!---------------------------------------------------------------------------
SUBROUTINE xc_( length, srd, svd, rho_in, ex_out, ec_out, vx_out, vc_out )
  !-------------------------------------------------------------------------
  !! Wrapper routine. Calls internal XC-driver routines or external ones
  !! from Libxc, depending on the input choice. See comments in routine 
  !! \(\texttt{xc}\) for variable explanations.
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_setting_params, ONLY: xc_func, xc_info, libxc_flags
#endif
  !
  USE kind_l,             ONLY: DP
  USE dft_setting_params, ONLY: iexch, icorr, is_libxc, rho_threshold_lda, &
                                finite_size_cell_volume_set
  USE qe_drivers_lda_lsda
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  INTEGER,  INTENT(IN) :: srd, svd
  REAL(DP), INTENT(IN) :: rho_in(length,srd)
  REAL(DP), INTENT(OUT) :: ex_out(length), vx_out(length,svd)
  REAL(DP), INTENT(OUT) :: ec_out(length), vc_out(length,svd)
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  INTEGER :: fkind_x, eflag
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
  REAL(DP), ALLOCATABLE :: zeta(:)
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
    ALLOCATE( rho_lxc(length*svd) )
    IF ( is_libxc(1) ) ALLOCATE( vx_lxc(length*svd) )
    IF ( is_libxc(2) ) ALLOCATE( vc_lxc(length*svd) )
    !
    ! ... set libxc input
    SELECT CASE( srd )
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
    IF (libxc_flags(1,0)==1) THEN
      CALL xc_f03_lda_exc_vxc( xc_func(1), lengthxc, rho_lxc(1), ex_out(1), vx_lxc(1) )
    ELSE
      CALL xc_f03_lda_vxc( xc_func(1), lengthxc, rho_lxc(1), vx_lxc(1) )
    ENDIF
  ENDIF
  !
  ! ... CORRELATION
  IF ( is_libxc(2) ) THEN
     CALL xc_f03_func_set_dens_threshold( xc_func(2), rho_threshold_lda )
     fkind_x = xc_f03_func_info_get_kind( xc_info(2) )
     IF (libxc_flags(2,0)==1) THEN
       CALL xc_f03_lda_exc_vxc( xc_func(2), lengthxc, rho_lxc(1), ec_out(1), vc_lxc(1) )
     ELSE
       CALL xc_f03_lda_vxc( xc_func(2), lengthxc, rho_lxc(1), vc_lxc(1) )
     ENDIF
  ENDIF
  !
  IF ( ((.NOT.is_libxc(1)) .OR. (.NOT.is_libxc(2))) &
        .AND. fkind_x/=XC_EXCHANGE_CORRELATION ) THEN
     !
     SELECT CASE( srd )
     CASE( 1 )
        !
        IF ((iexch==8 .AND. .NOT.is_libxc(1)) .OR. (icorr==10 .AND. &
            .NOT. is_libxc(2))) THEN
          IF (.NOT. finite_size_cell_volume_set) CALL xclib_error( 'XC',&
              'finite size corrected exchange used w/o initialization', 1 )
        ENDIF
        !$acc data copyin(rho_in), copyout(ex_out, ec_out, vx_out, vc_out)
        !$acc host_data use_device(rho_in, ex_out, ec_out, vx_out, vc_out)
        CALL xc_lda( length, rho_in(:,1), ex_out, ec_out, vx_out(:,1), vc_out(:,1) )
        !$acc end host_data
        !$acc end data
        !
     CASE( 2 )
        !
        ALLOCATE( zeta(length) )
        WHERE (ABS(rho_in(:,1)) > rho_threshold_lda) zeta(:) = rho_in(:,2) / ABS(rho_in(:,1))
        !$acc data copyin(rho_in, zeta), copyout(ex_out, ec_out, vx_out, vc_out)
        !$acc host_data use_device(rho_in, zeta, ex_out, ec_out, vx_out, vc_out)
        CALL xc_lsda( length, rho_in(:,1), zeta, ex_out, ec_out, vx_out, vc_out )
        !$acc end host_data
        !$acc end data
        DEALLOCATE( zeta )
        !
     CASE( 4 )
        !
        ALLOCATE( zeta(length) )
        WHERE (ABS(rho_in(:,1)) > rho_threshold_lda) zeta(:) = SQRT( rho_in(:,2)**2 + &
                                  rho_in(:,3)**2 + rho_in(:,4)**2 ) / ABS(rho_in(:,1)) ! amag/arho
        !$acc data copyin(rho_in, zeta), copyout(ex_out, ec_out, vx_out, vc_out)
        !$acc host_data use_device(rho_in, zeta, ex_out, ec_out, vx_out, vc_out)
        CALL xc_lsda( length, rho_in(:,1), zeta, ex_out, ec_out, vx_out, vc_out )
        !$acc end host_data
        !$acc end data
        DEALLOCATE( zeta )
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
  IF (svd == 1) THEN
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
     IF ( is_libxc(1) ) DEALLOCATE( vx_lxc )
     IF ( is_libxc(2) ) DEALLOCATE( vc_lxc )
  ENDIF
  !
#else
  !
  SELECT CASE( srd )
  CASE( 1 )
     !
     IF ((iexch==8 .AND. .NOT.is_libxc(1)) .OR. (icorr==10 .AND. &
            .NOT. is_libxc(2))) THEN
       IF (.NOT. finite_size_cell_volume_set) CALL xclib_error( 'XC',&
           'finite size corrected exchange used w/o initialization', 1 )
     ENDIF
     !
     !$acc data copyin(rho_in), copyout(ex_out, ec_out, vx_out, vc_out)
     !$acc host_data use_device(rho_in, ex_out, ec_out, vx_out, vc_out)
     CALL xc_lda( length, rho_in(:,1), ex_out, ec_out, vx_out(:,1), vc_out(:,1) )
     !$acc end host_data
     !$acc end data
     !
  CASE( 2 )
     !
     ALLOCATE( zeta(length) )
     WHERE (ABS(rho_in(:,1)) > rho_threshold_lda) zeta(:) = rho_in(:,2) / ABS(rho_in(:,1))
     !$acc data copyin(rho_in, zeta), copyout(ex_out, ec_out, vx_out, vc_out)
     !$acc host_data use_device(rho_in, zeta, ex_out, ec_out, vx_out, vc_out)
     CALL xc_lsda( length, rho_in(:,1), zeta, ex_out, ec_out, vx_out, vc_out )
     !$acc end host_data
     !$acc end data
     DEALLOCATE( zeta )
     ! 
  CASE( 4 )
     !
     ALLOCATE( zeta(length) )
     WHERE (ABS(rho_in(:,1)) > rho_threshold_lda) zeta(:) = SQRT( rho_in(:,2)**2 + &
                                rho_in(:,3)**2 + rho_in(:,4)**2 ) / ABS(rho_in(:,1)) ! amag/arho
     !$acc data copyin(rho_in, zeta), copyout(ex_out, ec_out, vx_out, vc_out)
     !$acc host_data use_device(rho_in, zeta, ex_out, ec_out, vx_out, vc_out)
     CALL xc_lsda( length, rho_in(:,1), zeta, ex_out, ec_out, vx_out, vc_out )
     !$acc end host_data
     !$acc end data
     DEALLOCATE( zeta )
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
END SUBROUTINE xc_
!
!---------------------------------------------------------------------------
SUBROUTINE xc_gpu( length, srd, svd, rho_in, ex_out, ec_out, vx_out, vc_out )
  !-------------------------------------------------------------------------
  !! Wrapper xc LDA - openACC version.  
  !! See comments in routine \(\texttt{xc}\) for variable explanations.
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_setting_params, ONLY: xc_func, xc_info, libxc_flags
#endif
  !
  USE kind_l,             ONLY: DP
  USE dft_setting_params, ONLY: iexch, icorr, is_libxc, rho_threshold_lda, &
                                finite_size_cell_volume_set
  USE qe_drivers_lda_lsda
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  INTEGER,  INTENT(IN) :: srd, svd
  REAL(DP), INTENT(IN) :: rho_in(length,srd)
  REAL(DP), INTENT(OUT) :: ex_out(length), vx_out(length,svd)
  REAL(DP), INTENT(OUT) :: ec_out(length), vc_out(length,svd)
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  LOGICAL :: is_libxc1, is_libxc2
  INTEGER :: fkind_x
  REAL(DP) :: amag
  REAL(DP), ALLOCATABLE :: rho_lxc(:)
  REAL(DP), ALLOCATABLE :: ex_lxc(:), ec_lxc(:)
  REAL(DP), ALLOCATABLE :: vx_lxc(:), vc_lxc(:)
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  REAL(DP), ALLOCATABLE :: zeta(:)
  REAL(DP) :: arho_ir
  INTEGER  :: ir
  !
  !$acc data deviceptr( rho_in(length,srd), ex_out(length), ec_out(length), &
  !$acc&                vx_out(length,svd), vc_out(length,svd) )
  !
#if defined(__LIBXC)
  !
  is_libxc1 = is_libxc(1)
  is_libxc2 = is_libxc(2)
  fkind_x = -1
  lengthxc = length
  !
  ALLOCATE( rho_lxc(length*svd) )
  ALLOCATE( ex_lxc(length), vx_lxc(length*svd) )
  ALLOCATE( ec_lxc(length), vc_lxc(length*svd) )
  !$acc data copyout( rho_lxc )
  !
  SELECT CASE( srd )  !... set libxc input
  CASE( 1 )
     !
     !$acc parallel loop
     DO ir = 1, length
       rho_lxc(ir) = ABS(rho_in(ir,1))
     ENDDO
     !
  CASE( 2 )
     !
     !$acc parallel loop
     DO ir = 1, length
        rho_lxc(2*ir-1) = (rho_in(ir,1) + rho_in(ir,2)) * 0.5_DP
        rho_lxc(2*ir)   = (rho_in(ir,1) - rho_in(ir,2)) * 0.5_DP
     ENDDO
     !
  CASE( 4 )
     !
     !$acc parallel loop
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
  !$acc end data
  !
  ! ... EXCHANGE
  IF ( is_libxc(1) ) THEN
    CALL xc_f03_func_set_dens_threshold( xc_func(1), rho_threshold_lda )
    IF (libxc_flags(1,0)==1) THEN
      CALL xc_f03_lda_exc_vxc( xc_func(1), lengthxc, rho_lxc(1), ex_out(1), vx_lxc(1) )
    ELSE
      CALL xc_f03_lda_vxc( xc_func(1), lengthxc, rho_lxc(1), vx_lxc(1) )
    ENDIF
  ENDIF
  !
  ! ... CORRELATION
  IF ( is_libxc(2) ) THEN
    CALL xc_f03_func_set_dens_threshold( xc_func(2), rho_threshold_lda )
    fkind_x = xc_f03_func_info_get_kind( xc_info(2) )
    IF (libxc_flags(2,0)==1) THEN
      CALL xc_f03_lda_exc_vxc( xc_func(2), lengthxc, rho_lxc(1), ec_out(1), vc_lxc(1) )
    ELSE
      CALL xc_f03_lda_vxc( xc_func(2), lengthxc, rho_lxc(1), vc_lxc(1) )
    ENDIF
  ENDIF
  !
  !
  IF ( ((.NOT.is_libxc(1)) .OR. (.NOT.is_libxc(2))) &
        .AND. fkind_x/=XC_EXCHANGE_CORRELATION ) THEN
     !
     SELECT CASE( srd )
     CASE( 1 )
        !
        IF (iexch==8 .OR. icorr==10) THEN
          IF (.NOT. finite_size_cell_volume_set) CALL xclib_error( 'XC',&
              'finite size corrected exchange used w/o initialization', 1 )
        ENDIF
        CALL xc_lda( length, rho_in(:,1), ex_out, ec_out, vx_out(:,1), vc_out(:,1) )
        !
     CASE( 2 )
        !
        ALLOCATE( zeta(length) )
        !$acc data create( zeta(length) )
        !$acc host_data use_device( zeta )
        !$acc parallel loop
        DO ir = 1, length
          arho_ir = ABS(rho_in(ir,1))
          IF (arho_ir > rho_threshold_lda) zeta(ir) = rho_in(ir,2) / arho_ir
        ENDDO
        CALL xc_lsda( length, rho_in(:,1), zeta, ex_out, ec_out, vx_out, vc_out )
        !$acc end host_data
        !$acc end data
        DEALLOCATE( zeta )
        !
     CASE( 4 )
        !
        ALLOCATE( zeta(length) )
        !$acc data create( zeta(length) )
        !$acc host_data use_device( zeta )
        !$acc parallel loop
        DO ir = 1, length
          arho_ir = ABS( rho_in(ir,1) )
          IF (arho_ir > rho_threshold_lda) zeta(ir) = SQRT( rho_in(ir,2)**2 + rho_in(ir,3)**2 + &
                                                          rho_in(ir,4)**2 ) / arho_ir ! amag/arho
        ENDDO
        CALL xc_lsda( length, rho_in(:,1), zeta, ex_out, ec_out, vx_out, vc_out )
        !$acc end host_data
        !$acc end data
        DEALLOCATE( zeta )
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
  !$acc data copyin( rho_lxc, ex_lxc, ec_lxc, vx_lxc, vc_lxc )
  !$acc parallel loop
  DO ir = 1, length
    IF (is_libxc1) ex_out(ir) = ex_lxc(ir)
    IF (is_libxc2) ec_out(ir) = ec_lxc(ir)
    IF (svd == 1) THEN
       IF (is_libxc1) vx_out(ir,1) = vx_lxc(ir)
       IF (is_libxc2) vc_out(ir,1) = vc_lxc(ir)
    ELSE
       IF (is_libxc1) THEN
         vx_out(ir,1) = vx_lxc(2*ir-1)
         vx_out(ir,2) = vx_lxc(2*ir)
       ENDIF
       IF (is_libxc2) THEN
         vc_out(ir,1) = vc_lxc(2*ir-1)
         vc_out(ir,2) = vc_lxc(2*ir)
       ENDIF
    ENDIF
  ENDDO
  !$acc end data
  DEALLOCATE( rho_lxc )
  DEALLOCATE( ex_lxc, vx_lxc )
  DEALLOCATE( ec_lxc, vc_lxc )
  !
#else
  !
  SELECT CASE( srd )
  CASE( 1 )
     !
     IF (iexch==8 .OR. icorr==10) THEN
       IF (.NOT. finite_size_cell_volume_set) CALL xclib_error( 'XC',&
           'finite size corrected exchange used w/o initialization', 1 )
     ENDIF
     !
     CALL xc_lda( length, rho_in(:,1), ex_out, ec_out, vx_out(:,1), vc_out(:,1) )
     !
  CASE( 2 )
     !
     ALLOCATE( zeta(length) )
     !$acc data create( zeta(length) )
     !$acc host_data use_device( zeta )
     !$acc parallel loop
     DO ir = 1, length
       arho_ir = ABS(rho_in(ir,1))
       IF (arho_ir > rho_threshold_lda) zeta(ir) = rho_in(ir,2) / arho_ir
     ENDDO
     CALL xc_lsda( length, rho_in(:,1), zeta, ex_out, ec_out, vx_out, vc_out )
     !$acc end host_data
     !$acc end data
     DEALLOCATE( zeta )
     !
   CASE( 4 )
     !
     ALLOCATE( zeta(length) )
     !$acc data create( zeta(length) )
     !$acc host_data use_device( zeta )
     !$acc parallel loop
     DO ir = 1, length
       arho_ir = ABS(rho_in(ir,1))
       IF (arho_ir > rho_threshold_lda) zeta(ir) = SQRT( rho_in(ir,2)**2 + rho_in(ir,3)**2 + &
                                                       rho_in(ir,4)**2 ) / arho_ir ! amag/arho
     ENDDO
     CALL xc_lsda( length, rho_in(:,1), zeta, ex_out, ec_out, vx_out, vc_out )
     !$acc end host_data
     !$acc end data
     DEALLOCATE( zeta )
     !
   CASE DEFAULT
     !
     CALL xclib_error( 'xc_LDA', 'Wrong ns input', 2 )
     !
  END SELECT
  !
#endif
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE xc_gpu
