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

  IMPLICIT NONE
  SAVE

  CHARACTER(LEN=256) :: fort_unit(30:44)
  ! ...  fort_unit = fortran units for saving physical quantity

CONTAINS


  SUBROUTINE printout_base_init( )

     USE io_files,  ONLY: tmp_dir, prefix
     USE io_global, ONLY: ionode, ionode_id
     USE mp_global, ONLY: intra_image_comm 
     USE mp, ONLY: mp_bcast

     INTEGER :: iunit, ierr, ios
     CHARACTER(LEN=256) :: pprefix
     ! ...  prefix combined with the output path

     pprefix = TRIM( tmp_dir ) // TRIM( prefix )

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
        fort_unit(41) = trim(pprefix)//'.spr'  ! wannier spread
        fort_unit(42) = trim(pprefix)//'.wfc'  ! wannier function
        fort_unit(43) = trim(pprefix)//'.hrs'  ! hirshfeld volumes 
        fort_unit(44) = trim(pprefix)//'.ncg'  ! number of cgsteps
        DO iunit = LBOUND( fort_unit, 1 ), UBOUND( fort_unit, 1 )
           OPEN(UNIT=iunit, FILE=fort_unit(iunit), &
               STATUS='unknown', POSITION='append', IOSTAT = ios )
           CLOSE( iunit )
           ierr = ierr + ABS(ios)
        END DO
     END IF

     CALL mp_bcast(ierr, ionode_id, intra_image_comm)
     IF( ierr /= 0 ) THEN
        CALL errore(' printout_base_init ', &
              ' error in opening files '//TRIM(pprefix)//'.XXX',ierr)
     END IF

    RETURN
  END SUBROUTINE printout_base_init


  SUBROUTINE printout_base_open( suffix )
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: suffix
    INTEGER :: iunit
    LOGICAL :: ok=.true.
    ! ...  Open units 30, 31, ... 44 for simulation output
    IF( PRESENT( suffix ) ) THEN
       IF( LEN( suffix ) /= 3 ) &
          CALL errore(" printout_base_open ", " wrong suffix ", 1 )
       ok = .false.
    END IF
    DO iunit = LBOUND( fort_unit, 1 ), UBOUND( fort_unit, 1 )
       IF( PRESENT( suffix ) ) THEN
          IF( index( fort_unit(iunit), suffix, back=.TRUE. ) == &
              ( len_trim( fort_unit(iunit) ) - 2 )                ) THEN
             OPEN( UNIT=iunit, FILE=fort_unit(iunit), STATUS='unknown', POSITION='append')
             ok = .true.
          END IF
       ELSE
          OPEN( UNIT=iunit, FILE=fort_unit(iunit), STATUS='unknown', POSITION='append')
       END IF
    END DO
    IF( PRESENT( suffix ) ) THEN
       IF( .NOT. ok ) &
          CALL errore(" printout_base_open ", " file with suffix "//suffix//" not found ", 1 )
    END IF
    RETURN
  END SUBROUTINE printout_base_open


  FUNCTION printout_base_unit( suffix )
    !   return the unit corresponding to a given suffix
    CHARACTER(LEN=*), INTENT(IN) :: suffix
    INTEGER :: printout_base_unit
    INTEGER :: iunit
    LOGICAL :: ok
    IF( LEN( suffix ) /= 3 ) &
       CALL errore(" printout_base_unit ", " wrong suffix ", 1 )
    ok = .false.
    DO iunit = LBOUND( fort_unit, 1 ), UBOUND( fort_unit, 1 )
       IF( index( fort_unit(iunit), suffix, back=.TRUE. ) == ( len_trim( fort_unit(iunit) ) - 2 ) ) THEN
          printout_base_unit = iunit
          ok = .true.
       END IF
    END DO
    IF( .NOT. ok ) &
       CALL errore(" printout_base_unit ", " file with suffix "//suffix//" not found ", 1 )
    RETURN
  END FUNCTION printout_base_unit


  FUNCTION printout_base_name( suffix )
    !  return the full name of a print out file with a given suffix
    CHARACTER(LEN=*), INTENT(IN) :: suffix
    CHARACTER(LEN=256) :: printout_base_name
    INTEGER :: iunit
    LOGICAL :: ok
    IF( LEN( suffix ) /= 3 ) &
       CALL errore(" printout_base_name ", " wrong suffix ", 1 )
    ok = .false.
    DO iunit = LBOUND( fort_unit, 1 ), UBOUND( fort_unit, 1 )
       IF( index( fort_unit(iunit), suffix, back=.TRUE. ) == ( len_trim( fort_unit(iunit) ) - 2 ) ) THEN
          printout_base_name = fort_unit(iunit)
          ok = .true.
       END IF
    END DO
    IF( .NOT. ok ) &
       CALL errore(" printout_base_name ", " file with suffix "//suffix//" not found ", 1 )
    RETURN
  END FUNCTION printout_base_name



  SUBROUTINE printout_base_close( suffix )
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: suffix
    INTEGER :: iunit
    LOGICAL :: topen
    LOGICAL :: ok
    ! ...   Close and flush unit 30, ... 44
    IF( PRESENT( suffix ) ) THEN
       IF( LEN( suffix ) /= 3 ) &
          CALL errore(" printout_base_close ", " wrong suffix ", 1 )
       ok = .false.
    END IF
    DO iunit = LBOUND( fort_unit, 1 ), UBOUND( fort_unit, 1 )
       IF( PRESENT( suffix ) ) THEN
          IF( index( fort_unit(iunit), suffix, back=.TRUE. ) == ( len_trim( fort_unit(iunit) ) - 2 ) ) THEN
             INQUIRE( UNIT=iunit, OPENED=topen )
             IF( topen ) CLOSE(iunit)
             ok = .true.
          END IF
       ELSE
          INQUIRE( UNIT=iunit, OPENED=topen )
          IF (topen) CLOSE(iunit)
       END IF
    END DO
    IF( PRESENT( suffix ) ) THEN
       IF( .NOT. ok ) &
          CALL errore(" printout_base_close ", " file with suffix "//suffix//" not found ", 1 )
    END IF
    RETURN
  END SUBROUTINE printout_base_close

  
  SUBROUTINE printout_pos( iunit, tau, nat, ityp, what, nfi, tps, label, fact, head )
    !
    USE kinds
    !
    INTEGER,          INTENT(IN)           :: iunit, nat, ityp(:)
    REAL(DP),        INTENT(IN)           :: tau( :, : )
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: what
    INTEGER,          INTENT(IN), OPTIONAL :: nfi
    REAL(DP),        INTENT(IN), OPTIONAL :: tps
    CHARACTER(LEN=3), INTENT(IN), OPTIONAL :: label( : )
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
    IF( PRESENT( label ) ) THEN
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
    !
END MODULE printout_base
