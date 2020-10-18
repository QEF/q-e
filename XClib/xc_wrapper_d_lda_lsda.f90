!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
SUBROUTINE dmxc( length, sr_d, rho_in, dmuxc )
  !---------------------------------------------------------------------
  !! Wrapper routine. Calls dmxc-driver routines from internal libraries
  !! or from the external one 'libxc', depending on the input choice.
  !
  ! Only two possibilities in the present version (LDA only):
  ! 1) iexch libxc + icorr libxc
  ! 2) iexch qe    + icorr qe
  !
  USE kind_l,            ONLY: DP
  USE dft_par_mod
  USE qe_drivers_d_lda_lsda
  !
#if defined(__LIBXC)
#include "xc_version.h"
  USE xc_f03_lib_m
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
  TYPE(xc_f03_func_t) :: xc_func
  TYPE(xc_f03_func_info_t) :: xc_info1, xc_info2
  INTEGER :: pol_unpol
  REAL(DP), ALLOCATABLE :: rho_lxc(:)
  REAL(DP), ALLOCATABLE :: dmxc_lxc(:), dmex_lxc(:), dmcr_lxc(:)
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
  IF ( (is_libxc(1) .OR. iexch==0) .AND. (is_libxc(2) .OR. icorr==0)) THEN
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
      CALL xclib_error( 'dmxc', 'The derivative of the xc potential with libxc &
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
    !
    ALLOCATE( dmex_lxc(length_dlxc), dmcr_lxc(length_dlxc), &
              dmxc_lxc(length_dlxc) )
    !
    ! ... DERIVATIVE FOR EXCHANGE
    dmex_lxc(:) = 0.0_DP
    IF (iexch /= 0) THEN    
       CALL xc_f03_func_init( xc_func, iexch, pol_unpol )
        xc_info1 = xc_f03_func_get_info( xc_func )
        CALL xc_f03_lda_fxc( xc_func, lengthxc, rho_lxc(1), dmex_lxc(1) )
       CALL xc_f03_func_end( xc_func )
    ENDIF    
    !
    ! ... DERIVATIVE FOR CORRELATION
    dmcr_lxc(:) = 0.0_DP
    IF (icorr /= 0) THEN
       CALL xc_f03_func_init( xc_func, icorr, pol_unpol )
        xc_info2 = xc_f03_func_get_info( xc_func )
        CALL xc_f03_lda_fxc( xc_func, lengthxc, rho_lxc(1), dmcr_lxc(1) )
       CALL xc_f03_func_end( xc_func )
    ENDIF
    !
    dmxc_lxc = (dmex_lxc + dmcr_lxc)*2.0_DP
    !
    IF (sr_d == 1) THEN
      dmuxc(:,1,1) = dmxc_lxc(:)
    ELSEIF (sr_d == 2) THEN
      DO ir = 1, length
        dmuxc(ir,1,1) = dmxc_lxc(3*ir-2)
        dmuxc(ir,1,2) = dmxc_lxc(3*ir-1)
        dmuxc(ir,2,1) = dmxc_lxc(3*ir-1)
        dmuxc(ir,2,2) = dmxc_lxc(3*ir)
      ENDDO
    ENDIF
    !
    DEALLOCATE( dmex_lxc, dmcr_lxc, dmxc_lxc )
    DEALLOCATE( rho_lxc )
    !
  ELSEIF ((.NOT.is_libxc(1)) .AND. (.NOT.is_libxc(2)) ) THEN
    !
    !CALL set_threshold_l( 'lda', small )
    rho_threshold_lda = small
    !
    IF ( sr_d == 1 ) CALL dmxc_lda( length, rho_in(:,1), dmuxc(:,1,1) )
    IF ( sr_d == 2 ) CALL dmxc_lsda( length, rho_in, dmuxc )
    !
  ELSE
    !
    CALL xclib_error( 'dmxc', 'Derivatives of exchange and correlation terms, &
                        & at present, must be both qe or both libxc.', 3 )
    !
  ENDIF
  !
#else
  !
  !CALL set_threshold_l( 'lda', small )
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
     CALL xclib_error( 'xc_LDA', 'Wrong ns input', 4 )
     !
  END SELECT
  !
#endif
  !
  !
  RETURN
  !
END SUBROUTINE dmxc
