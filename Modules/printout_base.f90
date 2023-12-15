!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! This module contains subroutines to print computed quantities to 
! standard output and ASCII file

MODULE printout_base

  USE kinds, only : DP
  IMPLICIT NONE
  SAVE

  INTEGER, EXTERNAL :: find_free_unit

CONTAINS


FUNCTION printout_base_name( suffix )
   !  return the full name of a print out file with a given suffix

   USE io_files,  ONLY: tmp_dir, prefix
   CHARACTER(LEN=*), INTENT(IN) :: suffix
   CHARACTER(LEN=256) :: printout_base_name
   printout_base_name = TRIM( tmp_dir ) // TRIM( prefix ) // TRIM(suffix)
   return
end function



  function printout_base_open( suffix )
    CHARACTER(LEN=*), INTENT(IN) :: suffix
    INTEGER :: printout_base_open
    LOGICAL :: ok=.true.

    printout_base_open = find_free_unit()
    if (printout_base_open <= 0 ) then
       call errore (" printout_base_open ", " cannot find a free unit ", 1)
    endif
    OPEN( UNIT=printout_base_open, FILE=printout_base_name(suffix),&
          STATUS='unknown', POSITION='append')
    RETURN
   END function printout_base_open

  SUBROUTINE printout_base_close( iunit )
    INTEGER, intent(in) :: iunit
    LOGICAL :: topen
    INQUIRE( UNIT=iunit, OPENED=topen )
    IF (topen) CLOSE(iunit)
  END SUBROUTINE printout_base_close

  
  SUBROUTINE printout_pos( iunit, tau, nat, ityp, what, nfi, tps, label, fact, head )
    !
    !
    INTEGER,          INTENT(IN)           :: iunit, nat
    INTEGER,          INTENT(IN), OPTIONAL :: ityp(:)
    REAL(DP),        INTENT(IN)           :: tau( :, : )
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: what
    INTEGER,          INTENT(IN), OPTIONAL :: nfi
    REAL(DP),        INTENT(IN), OPTIONAL :: tps
    CHARACTER(LEN=6), INTENT(IN), OPTIONAL :: label( : )
    REAL(DP),        INTENT(IN), OPTIONAL :: fact
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: head
    !
    INTEGER   :: ia, k
    REAL(DP) :: f
    !
    IF( PRESENT( fact ) ) THEN 
       f = fact
    ELSE
       f = 1.0_DP
    END IF
    !
    IF( PRESENT( head ) ) THEN
       WRITE( iunit, 10 ) head
    END IF
    !
    IF( PRESENT( what ) ) THEN
       IF ( what == 'xyz' ) WRITE( iunit, *) nat
    END IF
    !
    IF( PRESENT( nfi ) .AND. PRESENT( tps ) ) THEN
       WRITE( iunit, 30 ) nfi, tps
    ELSE IF( PRESENT( what ) ) THEN
       IF( what == 'pos' ) THEN
          WRITE( iunit, 40 )
       ELSE IF( what == 'vel' ) THEN
          WRITE( iunit, 50 )
       ELSE IF( what == 'for' ) THEN
          WRITE( iunit, 60 )
       END IF
    END IF
    !
    IF( PRESENT( label ) .and. PRESENT(ityp) ) THEN
       DO ia = 1, nat
         WRITE( iunit, 255 ) label(ityp(ia)), ( f * tau(k,ia),k = 1,3)
       END DO
    ELSE
       DO ia = 1, nat
         WRITE( iunit, 252 ) (tau(k,ia),k = 1,3)
       END DO
    END IF
 10 FORMAT(3X,A)
 30 FORMAT(I8,1X,F13.8)
 40 FORMAT(3X,'ATOMIC_POSITIONS')
 50 FORMAT(3X,'ATOMIC_VELOCITIES')
 60 FORMAT(3X,'Forces acting on atoms (au):')
255 FORMAT(3X,A3,3E25.14)
252 FORMAT(3E25.14)
    RETURN
  END SUBROUTINE printout_pos

 

  SUBROUTINE printout_cell( iunit, h, nfi, tps )
    !
    USE kinds
    !
    INTEGER,   INTENT(IN)           :: iunit
    REAL(DP), INTENT(IN)           :: h(3,3)
    INTEGER,   INTENT(IN), OPTIONAL :: nfi
    REAL(DP), INTENT(IN), OPTIONAL :: tps
    !
    INTEGER :: i, j
    !
    IF( PRESENT( nfi ) .AND. PRESENT( tps ) ) THEN
       WRITE( iunit, 30 ) nfi, tps
    ELSE
       WRITE( iunit, 40 )
    END IF
    !
    DO i = 1, 3
       WRITE( iunit, 100 ) (h(i,j),j=1,3)
    END DO
    !
 30 FORMAT(I8,1X,F13.8)
 40 FORMAT(3X,'CELL_PARAMETERS')
100 FORMAT(3F14.8)
    RETURN
  END SUBROUTINE printout_cell



  SUBROUTINE printout_stress( iunit, str, nfi, tps )
    !
    USE kinds
    !
    INTEGER,   INTENT(IN)           :: iunit
    REAL(DP), INTENT(IN)           :: str(3,3)
    INTEGER,   INTENT(IN), OPTIONAL :: nfi
    REAL(DP), INTENT(IN), OPTIONAL :: tps
    !
    INTEGER :: i, j
    !
    IF( PRESENT( nfi ) .AND. PRESENT( tps ) ) THEN
       WRITE( iunit, 30 ) nfi, tps
    ELSE
       WRITE( iunit, 40 )
    END IF
    !
    DO i = 1, 3
       WRITE( iunit, 100 ) (str(i,j),j=1,3)
    END DO
    !
 30 FORMAT(I8,1X,F13.8)
 40 FORMAT(3X,'Total stress (GPa)')
100 FORMAT(3(F18.8,1X))
    RETURN
  END SUBROUTINE printout_stress

  SUBROUTINE printout_vefftsvdw( iunit, veff, nat, nfi, tps )
    !
    USE kinds
    !
    INTEGER,   INTENT(IN)           :: iunit, nat
    REAL(DP), INTENT(IN)           :: veff(nat)
    INTEGER,   INTENT(IN), OPTIONAL :: nfi
    REAL(DP), INTENT(IN), OPTIONAL :: tps
    !
    INTEGER :: i, j
    !
    IF( PRESENT( nfi ) .AND. PRESENT( tps ) ) THEN
       WRITE( iunit, 30 ) nfi, tps
    ELSE
       WRITE( iunit, 40 )
    END IF
    !
    DO i = 1, nat 
       WRITE( iunit, 100 ) veff(i) 
    END DO
    !
 30 FORMAT(I8,1X,F13.8)
 40 FORMAT(3X,'Veff tsvdw')
100 FORMAT(F20.10)
    RETURN
  END SUBROUTINE printout_vefftsvdw

  SUBROUTINE printout_wfc( iunit, wfc_temp, nband, nfi, tps, iss )
    !
    USE kinds
    !
    INTEGER,   INTENT(IN)           :: iunit, nband 
    REAL(DP), INTENT(IN)           :: wfc_temp(3,nband)
    INTEGER,   INTENT(IN)           :: nfi
    REAL(DP), INTENT(IN)           :: tps
    INTEGER, INTENT(IN), OPTIONAL  :: iss 
    !
    INTEGER :: i, j
    !
    IF( PRESENT( iss ) ) THEN
       WRITE( iunit, 40 ) nfi, tps, iss
    ELSE
       WRITE( iunit, 30 ) nfi, tps
    END IF
    !
    DO i = 1, nband 
       WRITE( iunit, 100 ) (wfc_temp(j,i),j=1,3) 
    END DO
    !
 30 FORMAT(I8,1X,F13.8)
 40 FORMAT(I7,1X,F11.8,1X,"spin=",I5)
100 FORMAT(3E25.14)
    RETURN
  END SUBROUTINE printout_wfc
    !------------------------------------------------------------------------
    SUBROUTINE save_print_counter( iter, wunit )
      !------------------------------------------------------------------------
      !
      ! ... a counter indicating the last successful printout iteration is saved
      !
      USE io_global, ONLY: ionode, ionode_id
      USE io_files, ONLY : iunpun, create_directory, restart_dir
      USE mp, ONLY: mp_bcast
      USE mp_images, ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(IN) :: iter
      INTEGER,          INTENT(IN) :: wunit
      !
      INTEGER            :: ierr
      CHARACTER(LEN=256) :: filename, dirname
      !
      !
      dirname = restart_dir( wunit )
      !
      CALL create_directory( TRIM( dirname ) )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // 'print_counter'
         !
         OPEN( UNIT = iunpun, FILE = filename, FORM = 'formatted', &
                 STATUS = 'unknown', IOSTAT = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'save_print_counter', &
                   'cannot open restart file for writing', ierr )
      !
      IF ( ionode ) THEN
         !
         WRITE ( iunpun, '("LAST SUCCESSFUL PRINTOUT AT STEP:",/,i5 )' ) iter
         !
         CLOSE ( iunpun, STATUS = 'keep' )
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE save_print_counter
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_print_counter( nprint_nfi, runit )
      !------------------------------------------------------------------------
      !
      ! ... the counter indicating the last successful printout iteration 
      ! ... is read here
      !
      USE io_global, ONLY: ionode, ionode_id
      USE io_files, ONLY : iunpun, restart_dir
      USE mp, ONLY: mp_bcast
      USE mp_images, ONLY : intra_image_comm
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: nprint_nfi
      INTEGER,          INTENT(IN)  :: runit
      !
      INTEGER            :: ierr
      CHARACTER(LEN=256) :: filename, dirname
      !
      !
      dirname = restart_dir( runit )
      !
      IF ( ionode ) THEN
         !
         filename = TRIM( dirname ) // 'print_counter'
         !
         OPEN( UNIT = iunpun, FILE = filename, FORM = 'formatted', &
                 STATUS = 'old', IOSTAT = ierr )
         !
         IF ( ierr > 0 ) THEN
            !
            nprint_nfi = -1
            !
         ELSE
            !
            READ ( iunpun, * )
            READ ( iunpun, * ) nprint_nfi
            !
            CLOSE ( iunpun, STATUS = 'keep' )
            !
         END IF
         !
      END IF
      !
      CALL mp_bcast( nprint_nfi, ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_print_counter   


   SUBROUTINE print_eigenvalues( ei_unit, tfile, tstdout, nfi, tps, nspin, ei, nupdwn )
      !
      use constants,  only : autoev 
      USE io_global,  ONLY : stdout, ionode
      !
      INTEGER,  INTENT(IN) :: ei_unit, nupdwn(:)
      LOGICAL,  INTENT(IN) :: tfile, tstdout
      INTEGER,  INTENT(IN) :: nfi, nspin
      REAL(DP), INTENT(IN) :: tps, ei(:,:)
      !
      INTEGER :: i, j, ik
      !
      IF ( tfile ) THEN
          WRITE(ei_unit,30) nfi, tps
      END IF
      !
      ik = 1
      !
      DO j = 1, nspin
         !
         IF( tstdout ) THEN
            WRITE( stdout,1002) ik, j
            WRITE( stdout,1004) ( ei( i, j ) * autoev, i = 1, nupdwn(j) )
         END IF
         !
         IF( tfile ) THEN
            WRITE(ei_unit,1010) ik, j
            WRITE(ei_unit,1020) ( ei( i, j ) * autoev, i = 1, nupdwn(j) )
         END IF
         !
      END DO
      !
  30  FORMAT(2X,'STEP:',I7,1X,F10.2)
 1002 FORMAT(/,3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2,/)
 1004 FORMAT(10F8.2)
 1010 FORMAT(3X,'Eigenvalues (eV), kp = ',I3, ' , spin = ',I2)
 1020 FORMAT(10F8.2)
      !
      RETURN
   END SUBROUTINE print_eigenvalues

    !
END MODULE printout_base
