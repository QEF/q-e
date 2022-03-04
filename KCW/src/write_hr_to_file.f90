!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
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
  USE control_kcw,          ONLY : do_bands, Hamlt_R, num_wann, irvect,&
                                   num_wann_occ, num_wann_emp
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
  INTEGER :: ir, i
  INTEGER :: iband, jband
  INTEGER :: iband_, jband_
  INTEGER :: nwann
  INTEGER :: ifile, fileunit
  CHARACTER(LEN=256) :: filename
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
  DO ifile = 1, 3
    !
    CALL date_and_tim( cdate, ctime )
    header = 'Written on '//cdate//' at '//ctime
    !
    fileunit = 99 + ifile
    !
    IF ( fileunit == 100 ) THEN
      !
      WRITE( filename, 20 ) TRIM(prefix)
      nwann = num_wann
      !
    ELSE IF ( fileunit == 101 ) THEN
      !
      WRITE( filename, 21 ) TRIM(prefix), 'occ'
      nwann = num_wann_occ
      !
    ELSE
      !
      WRITE( filename, 21 ) TRIM(prefix), 'emp'
      nwann = num_wann_emp
      !
    ENDIF
    !
    IF ( ionode ) THEN
      !
      OPEN( fileunit, file=TRIM(filename), status='unknown' )
      !
      WRITE( fileunit, * ) header
      !
      IF ( fileunit == 100 ) THEN
        !
        WRITE( fileunit, '(4x,"Rtot =",i5,6x,"num_wann =",i5)' ) nkstot/nspin, nwann
        !
      ELSE
        !
        WRITE( fileunit, * ) nwann
        WRITE( fileunit, * ) nkstot/nspin
        WRITE( fileunit, '(15I5)' ) (1, i=1, nkstot/nspin)
        !
      ENDIF
      !
      DO ir = 1, nkstot/nspin
        DO iband = 1, nwann
          DO jband = 1, nwann
            !
            IF ( fileunit == 102 ) THEN
              iband_ = iband + num_wann_occ
              jband_ = jband + num_wann_occ
            ELSE
              iband_ = iband
              jband_ = jband
            ENDIF
            !
            WRITE( fileunit, '(5i5,2f12.6)' ) irvect(:, ir), jband, iband, Hamlt_R(ir, jband_, iband_) * rytoev
            !
          ENDDO
        ENDDO
      ENDDO
      !
      CLOSE( fileunit )
      !
    ENDIF
    !
  ENDDO
  !
20 FORMAT( A, '.kcw_hr.dat' )
21 FORMAT( A, '.kcw_hr_', A, '.dat' )
  !
  !
END SUBROUTINE write_hr_to_file
