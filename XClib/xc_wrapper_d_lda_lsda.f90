!
! Copyright (C) 2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------
SUBROUTINE dmxc( length, sr_d, rho_in, dmuxc )
  !----------------------------------------------------------------------
  !! Wrapper routine. Calls internal dmxc-driver routines or the external
  !! ones from Libxc, depending on the input choice.
  !
  USE kind_l,          ONLY: DP
  USE dft_par_mod,     ONLY: iexch, icorr, is_libxc, rho_threshold_lda
  USE qe_drivers_d_lda_lsda
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
  USE dft_par_mod,     ONLY: xc_func, xc_info
#endif
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: length
  !! length of the I/O arrays
  INTEGER,  INTENT(IN) :: sr_d
  !! number of spin components
  REAL(DP), INTENT(IN) :: rho_in(length,sr_d)
  !! charge density
  REAL(DP), INTENT(OUT) :: dmuxc(length,sr_d,sr_d)
  !! the derivative of the xc potential
  !
  ! ... local variables
  !
#if defined(__LIBXC)
  INTEGER :: pol_unpol, fkind_x
  REAL(DP), ALLOCATABLE :: rho_lxc(:)
  REAL(DP), ALLOCATABLE :: dmex_lxc(:), dmcr_lxc(:)
  LOGICAL :: exch_lxc_avail, corr_lxc_avail
#if (XC_MAJOR_VERSION > 4)
  INTEGER(8) :: lengthxc
#else
  INTEGER :: lengthxc
#endif
#endif
  !
  INTEGER :: ir, length_lxc, length_dlxc
  REAL(DP), PARAMETER :: small = 1.E-10_DP, rho_trash = 0.5_DP
  !
#if defined(__LIBXC)
  !
  lengthxc = length
  !
  IF ( ANY(is_libxc(1:2)) ) THEN
    !
    length_lxc = length*sr_d
    !
    ! ... set libxc input
    SELECT CASE( sr_d )
    CASE( 1 )
      !
      ALLOCATE( rho_lxc(length_lxc) )
      pol_unpol = 1
      rho_lxc = rho_in(:,1) 
      !
    CASE( 2 )
      !
      ALLOCATE( rho_lxc(length_lxc) )
      pol_unpol = 2
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
  ENDIF
  !
  IF ( is_libxc(1) ) THEN
    ALLOCATE( dmex_lxc(length_dlxc) )
    ! ... DERIVATIVE FOR EXCHANGE
    dmex_lxc(:) = 0.0_DP
    IF (iexch /= 0) THEN    
       CALL xc_f03_lda_fxc( xc_func(1), lengthxc, rho_lxc(1), dmex_lxc(1) )
    ENDIF  
  ENDIF
  !
  IF ( is_libxc(2) ) THEN
    ALLOCATE( dmcr_lxc(length_dlxc) )
    ! ... DERIVATIVE FOR CORRELATION
    dmcr_lxc(:) = 0.0_DP
    IF (icorr /= 0) THEN
       fkind_x  = xc_f03_func_info_get_kind( xc_info(2) )
       CALL xc_f03_lda_fxc( xc_func(2), lengthxc, rho_lxc(1), dmcr_lxc(1) )
    ENDIF
  ENDIF
  !
  dmuxc = 0.0_DP
  !
  IF ( ((.NOT.is_libxc(1)) .OR. (.NOT.is_libxc(2))) &
        .AND. fkind_x/=XC_EXCHANGE_CORRELATION ) THEN
    !
    rho_threshold_lda = small
    !
    IF ( sr_d == 1 ) CALL dmxc_lda( length, rho_in(:,1), dmuxc(:,1,1) )
    IF ( sr_d == 2 ) CALL dmxc_lsda( length, rho_in, dmuxc )
    IF ( sr_d == 4 ) CALL dmxc_nc( length, rho_in(:,1), rho_in(:,2:4), dmuxc )
    !
  ENDIF
  !
  !
  IF ( ANY(is_libxc(1:2)) ) THEN
    SELECT CASE( sr_d )
    CASE( 1 )
      !
      IF ( is_libxc(1) ) dmuxc(:,1,1) = dmuxc(:,1,1) + dmex_lxc(:)*2.0_DP
      IF ( is_libxc(2) ) dmuxc(:,1,1) = dmuxc(:,1,1) + dmcr_lxc(:)*2.0_DP
      !
    CASE( 2 )
      !
      IF ( is_libxc(1) ) THEN
        DO ir = 1, length
          dmuxc(ir,1,1) = dmuxc(ir,1,1) + dmex_lxc(3*ir-2)*2.0_DP
          dmuxc(ir,1,2) = dmuxc(ir,1,2) + dmex_lxc(3*ir-1)*2.0_DP
          dmuxc(ir,2,1) = dmuxc(ir,2,1) + dmex_lxc(3*ir-1)*2.0_DP
          dmuxc(ir,2,2) = dmuxc(ir,2,2) + dmex_lxc(3*ir)  *2.0_DP
        ENDDO
        DEALLOCATE( dmex_lxc )
      ENDIF
      !
      IF ( is_libxc(2) ) THEN
        DO ir = 1, length  
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
#else
  !
  rho_threshold_lda = small
  !
  SELECT CASE( sr_d )
  CASE( 1 )
     !
     CALL dmxc_lda( length, rho_in(:,1), dmuxc(:,1,1) )
     !
  CASE( 2 )
     !
     CALL dmxc_lsda( length, rho_in, dmuxc )
     ! 
  CASE( 4 )
     !
     CALL dmxc_nc( length, rho_in(:,1), rho_in(:,2:4), dmuxc )
     !
  CASE DEFAULT
     !
     CALL xclib_error( 'dmxc', 'Wrong ns input', 4 )
     !
  END SELECT
  !
#endif
  !
  !
  RETURN
  !
END SUBROUTINE dmxc
