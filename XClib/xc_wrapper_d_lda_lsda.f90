!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE dmxc( length, srd, rho_in, dmuxc, gpu_args_ )
  !----------------------------------------------------------------------
  !! Wrapper routine. Calls internal dmxc-driver routines or the external
  !! ones from Libxc, depending on the input choice.
  !
  USE kind_l,   ONLY: DP
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: srd
  !! number of spin components
  REAL(DP), INTENT(IN) :: rho_in(length,srd)
  !! charge density
  REAL(DP), INTENT(OUT) :: dmuxc(length,srd,srd)
  !! the derivative of the xc potential
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
    !$acc data present( rho_in, dmuxc )
    CALL dmxc_( length, srd, rho_in, dmuxc )
    !$acc end data
    !
  ELSE
    !
    !$acc data copyin( rho_in ), copyout( dmuxc )
    CALL dmxc_( length, srd, rho_in, dmuxc )
    !$acc end data
    !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE
!
!------------------------------------------------------------------------
SUBROUTINE dmxc_( length, srd, rho_in, dmuxc )
  !----------------------------------------------------------------------
  !! Wrapper routine. Calls internal dmxc-driver routines or the external
  !! ones from Libxc, depending on the input choice.
  !
  USE kind_l,               ONLY: DP
  USE dft_setting_params,   ONLY: iexch, icorr, is_libxc, rho_threshold_lda
  USE qe_drivers_d_lda_lsda
  !
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
  INTEGER,  INTENT(IN) :: srd
  !! number of spin components
  REAL(DP), INTENT(IN) :: rho_in(length,srd)
  !! charge density
  REAL(DP), INTENT(OUT) :: dmuxc(length,srd,srd)
  !! the derivative of the xc potential
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  INTEGER :: pol_unpol
  REAL(DP), ALLOCATABLE :: rho_lxc(:)
  REAL(DP), ALLOCATABLE :: dmex_lxc(:), dmcr_lxc(:)
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  !
  LOGICAL :: fkind_is_XC
  INTEGER :: ir, length_lxc, length_dlxc, fkind_x
  REAL(DP), PARAMETER :: small=1.E-10_DP, rho_trash=0.5_DP
  !
  !$acc data present( rho_in, dmuxc )
  !
  !$acc kernels
  dmuxc(:,:,:) = 0.0_DP
  !$acc end kernels
  !
  fkind_x = -1
  fkind_is_XC = .FALSE.
  !
#if defined(__LIBXC)
  !
  lengthxc = length
  !
  IF ( ANY(is_libxc(1:2)) ) THEN
    !
    length_lxc = length*srd
    !
    ALLOCATE( rho_lxc(length_lxc) )
    !$acc data copyout( rho_lxc )
    !
    ! ... set libxc input
    SELECT CASE( srd )
    CASE( 1 )
      !
      pol_unpol = 1
      !$acc parallel loop
      DO ir = 1, length
        rho_lxc(ir) = rho_in(ir,1)
      ENDDO
      !
    CASE( 2 )
      !
      pol_unpol = 2
      !$acc parallel loop
      DO ir = 1, length
        rho_lxc(2*ir-1) = rho_in(ir,1)
        rho_lxc(2*ir)   = rho_in(ir,2)
      ENDDO
      !
    CASE( 4 )
      !
      CALL xclib_error( 'dmxc', 'The Libxc derivative of the XC potential &
                           &is not available for noncollinear case', 1 )
      !
    CASE DEFAULT
      !
      CALL xclib_error( 'dmxc', 'Wrong number of spin dimensions', 2 )
      !
    END SELECT
    !
    length_dlxc = length
    IF (pol_unpol == 2) length_dlxc = length*3
    !
    !$acc end data
    !
  ENDIF
  !
  IF ( is_libxc(1) ) THEN
    ALLOCATE( dmex_lxc(length_dlxc) )
    ! ... DERIVATIVE FOR EXCHANGE
    dmex_lxc(:) = 0.0_DP
    IF (iexch /= 0) THEN
       CALL xc_f03_func_set_dens_threshold( xc_func(1), rho_threshold_lda )
       CALL xc_f03_lda_fxc( xc_func(1), lengthxc, rho_lxc(1), dmex_lxc(1) )
    ENDIF  
  ENDIF
  !
  fkind_x = -100
  IF ( is_libxc(2) ) THEN
    ALLOCATE( dmcr_lxc(length_dlxc) )
    ! ... DERIVATIVE FOR CORRELATION
    dmcr_lxc(:) = 0.0_DP
    IF (icorr /= 0) THEN
       fkind_x  = xc_f03_func_info_get_kind( xc_info(2) )
       CALL xc_f03_func_set_dens_threshold( xc_func(2), rho_threshold_lda )
       CALL xc_f03_lda_fxc( xc_func(2), lengthxc, rho_lxc(1), dmcr_lxc(1) )
    ENDIF
  ENDIF
  !
  fkind_is_XC = (fkind_x==XC_EXCHANGE_CORRELATION)
  !
#endif
  !
  IF ( ((.NOT.is_libxc(1)) .OR. (.NOT.is_libxc(2))) &
        .AND. (.NOT.fkind_is_XC) ) THEN
    rho_threshold_lda = small
    IF ( srd == 1 ) CALL dmxc_lda( length, rho_in(:,1), dmuxc(:,1,1) )
    IF ( srd == 2 ) CALL dmxc_lsda( length, rho_in, dmuxc )
    IF ( srd == 4 ) CALL dmxc_nc( length, rho_in, dmuxc )
  ENDIF
  !
#if defined(__LIBXC)
  !
  IF ( ANY(is_libxc(1:2)) ) THEN
    SELECT CASE( srd )
    CASE( 1 )
      !
      IF ( is_libxc(1) ) THEN
        !$acc parallel loop copyin( dmex_lxc )
        DO ir = 1, length
          IF (rho_in(ir,1)<=rho_threshold_lda ) CYCLE
          dmuxc(ir,1,1) = dmuxc(ir,1,1) + dmex_lxc(ir)*2.0_DP
        ENDDO
        DEALLOCATE( dmex_lxc )
      ENDIF
      !
      IF ( is_libxc(2) ) THEN
        !$acc parallel loop copyin( dmcr_lxc )
        DO ir = 1, length
          IF (rho_in(ir,1)<=rho_threshold_lda ) CYCLE
          dmuxc(ir,1,1) = dmuxc(ir,1,1) + dmcr_lxc(ir)*2.0_DP
        ENDDO
        DEALLOCATE( dmcr_lxc )
      ENDIF
      !
    CASE( 2 )
      !
      IF ( is_libxc(1) ) THEN
        !$acc parallel loop copyin( dmex_lxc )
        DO ir = 1, length
          IF (rho_in(ir,1)<=rho_threshold_lda .OR. &
              rho_in(ir,2)<=rho_threshold_lda) CYCLE
          dmuxc(ir,1,1) = dmuxc(ir,1,1) + dmex_lxc(3*ir-2)*2.0_DP
          dmuxc(ir,1,2) = dmuxc(ir,1,2) + dmex_lxc(3*ir-1)*2.0_DP
          dmuxc(ir,2,1) = dmuxc(ir,2,1) + dmex_lxc(3*ir-1)*2.0_DP
          dmuxc(ir,2,2) = dmuxc(ir,2,2) + dmex_lxc(3*ir)  *2.0_DP
        ENDDO
        DEALLOCATE( dmex_lxc )
      ENDIF
      !
      IF ( is_libxc(2) ) THEN
        !$acc parallel loop copyin( dmcr_lxc )
        DO ir = 1, length
          IF (rho_in(ir,1)<=rho_threshold_lda .OR. &
              rho_in(ir,2)<=rho_threshold_lda) CYCLE
          dmuxc(ir,1,1) = dmuxc(ir,1,1) + dmcr_lxc(3*ir-2)*2.0_DP
          dmuxc(ir,1,2) = dmuxc(ir,1,2) + dmcr_lxc(3*ir-1)*2.0_DP
          dmuxc(ir,2,1) = dmuxc(ir,2,1) + dmcr_lxc(3*ir-1)*2.0_DP
          dmuxc(ir,2,2) = dmuxc(ir,2,2) + dmcr_lxc(3*ir)  *2.0_DP
        ENDDO
        DEALLOCATE( dmcr_lxc )
      ENDIF
      !
    END SELECT
  ENDIF
  !
  IF ( ANY(is_libxc(1:2)) ) DEALLOCATE( rho_lxc )
  !
#endif
  !
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE dmxc_
