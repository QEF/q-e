!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! This module contains subroutines to print computed quantities to 
! standard output and ASCII file

MODULE printout_base

  IMPLICIT NONE
  SAVE

  CHARACTER(LEN=75) :: title
  ! ...  title of the simulation

  CHARACTER(LEN=256) :: fort_unit(30:40)
  ! ...  fort_unit = fortran units for saving physical quantity

  CHARACTER(LEN=256) :: pprefix
  ! ...  prefix combined with the output path

CONTAINS


  SUBROUTINE printout_base_init( outdir, prefix )

     USE io_global, ONLY: ionode, ionode_id
     USE mp_global, ONLY: group
     USE mp, ONLY: mp_bcast

     INTEGER :: iunit, ierr
     CHARACTER(LEN=*), INTENT(IN) :: outdir
     CHARACTER(LEN=*), INTENT(IN) :: prefix
     CHARACTER(LEN=256) :: file_name


     IF( prefix /= ' ' ) THEN
        pprefix = TRIM( prefix )
     ELSE
        pprefix = 'fpmd'
     END IF

     IF( outdir /= ' ' ) THEN
        pprefix = TRIM( outdir ) // '/' // TRIM( pprefix )
     END IF

     ierr = 0

     IF( ionode ) THEN
        fort_unit(30) = trim(pprefix)//'.con'
        fort_unit(31) = trim(pprefix)//'.eig'
        fort_unit(32) = trim(pprefix)//'.pol'
        fort_unit(33) = trim(pprefix)//'.evp'
        fort_unit(34) = trim(pprefix)//'.vel'
        fort_unit(35) = trim(pprefix)//'.pos'
        fort_unit(36) = trim(pprefix)//'.cel'
        fort_unit(37) = trim(pprefix)//'.for'
        fort_unit(38) = trim(pprefix)//'.str'
        fort_unit(39) = trim(pprefix)//'.nos'
        fort_unit(40) = trim(pprefix)//'.the'
        DO iunit = LBOUND( fort_unit, 1 ), UBOUND( fort_unit, 1 )
           OPEN(UNIT=iunit, FILE=fort_unit(iunit), &
               STATUS='unknown', POSITION='append', IOSTAT = ierr )
           CLOSE( iunit )
        END DO
     END IF

     CALL mp_bcast(ierr, ionode_id, group)
     IF( ierr /= 0 ) &
        CALL errore(' printout_base_init ',' error in opening unit, check outdir ',iunit)

    RETURN
  END SUBROUTINE printout_base_init


  SUBROUTINE printout_base_open( )
    INTEGER :: iunit
    ! ...  Open units 30, 31, ... 40 for simulation output
    DO iunit = LBOUND( fort_unit, 1 ), UBOUND( fort_unit, 1 )
       OPEN( UNIT=iunit, FILE=fort_unit(iunit), STATUS='unknown', POSITION='append')
    END DO
    RETURN
  END SUBROUTINE printout_base_open

  SUBROUTINE printout_base_close( )
    INTEGER :: iunit
    LOGICAL :: topen
    ! ...   Close and flush unit 30, ... 40
    DO iunit = LBOUND( fort_unit, 1 ), UBOUND( fort_unit, 1 )
       INQUIRE( UNIT=iunit, OPENED=topen )
       IF (topen) THEN
          CLOSE(iunit)
       END IF
    END DO
    RETURN
  END SUBROUTINE printout_base_close

  
  SUBROUTINE printout_pos( iunit, tau, nat, what, nfi, tps, label, fact, sort )
    !
    USE kinds
    !
    INTEGER,          INTENT(IN)           :: iunit, nat
    REAL(DP),        INTENT(IN)           :: tau( :, : )
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: what
    INTEGER,          INTENT(IN), OPTIONAL :: nfi
    REAL(DP),        INTENT(IN), OPTIONAL :: tps
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: label( : )
    REAL(DP),        INTENT(IN), OPTIONAL :: fact
    INTEGER,          INTENT(IN), OPTIONAL :: sort( : )
    !
    INTEGER   :: ia, k
    REAL(DP) :: f
    !
    IF( PRESENT( fact ) ) THEN 
       f = fact
    ELSE
       f = 1.0d0
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
    IF( PRESENT( label ) ) THEN
       IF( PRESENT( sort ) ) THEN
         DO ia = 1, nat
           WRITE( iunit, 255 ) label( sort(ia) ), ( f * tau(k, sort(ia) ),k = 1,3)
         END DO
       ELSE
         DO ia = 1, nat
           WRITE( iunit, 255 ) label(ia), ( f * tau(k,ia),k = 1,3)
         END DO
       END IF
    ELSE
       DO ia = 1, nat
         WRITE( iunit, 252 ) (tau(k,ia),k = 1,3)
       END DO
    END IF
 30 FORMAT(I7,1X,F11.8)
 40 FORMAT(3X,'ATOMIC_POSITIONS')
 50 FORMAT(3X,'ATOMIC_VELOCITIES')
 60 FORMAT(3X,'Forces acting on atoms (au):')
255 FORMAT(3X,A3,3E14.6)
252 FORMAT(3E14.6)
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
 30 FORMAT(I7,1X,F11.8)
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
 30 FORMAT(I7,1X,F11.8)
 40 FORMAT(3X,'Total stress (GPa)')
100 FORMAT(3(F18.8,1X))
    RETURN
  END SUBROUTINE printout_stress

END MODULE printout_base
