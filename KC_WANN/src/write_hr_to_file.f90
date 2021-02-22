!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_hr_to_file()
  !---------------------------------------------------------------------
  !
  !! This routine calculates KC matrix H(R) - if not calculated
  !! already - and it prints it into a formatted file
  !
  USE control_kc_wann,      ONLY : do_bands, Hamlt_R, num_wann, irvect
  USE klist,                ONLY : nkstot
  USE lsda_mod,             ONLY : nspin
  USE interpolation,        ONLY : real_ham
  USE io_files,             ONLY : prefix
  USE io_global,            ONLY : ionode
  USE constants,            ONLY : rytoev
  !
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=9)  :: cdate, ctime
  CHARACTER(LEN=33) :: header
  INTEGER :: ir
  INTEGER :: iband, jband
  !
  !
  IF ( .NOT. do_bands ) THEN
    !
    ALLOCATE( Hamlt_R(nkstot/nspin,num_wann,num_wann) )
    !
    CALL real_ham( Hamlt_R )
    !
  ENDIF
  !
  CALL date_and_tim( cdate, ctime )
  header = 'Written on '//cdate//' at '//ctime
  !
  !
  IF ( ionode ) THEN
    !
    OPEN( 100, file=TRIM(prefix)//'.kc_hr.dat', status='unknown' )
    !
    WRITE( 100, * ) header
    WRITE( 100, '(4x,"Rtot =",i5,6x,"num_wann =",i5)' ) nkstot/nspin, num_wann
    !
    DO ir = 1, nkstot/nspin
      DO iband = 1, num_wann
        DO jband = 1, num_wann
          !
          WRITE( 100, '(5i5,2f12.6)' ) irvect(:,ir), jband, iband, Hamlt_R(ir,jband,iband)*rytoev
          !
        ENDDO
      ENDDO
    ENDDO
    !
    CLOSE( 100 )
    !
  ENDIF
  !
  !
END SUBROUTINE write_hr_to_file
