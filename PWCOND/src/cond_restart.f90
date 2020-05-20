!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE cond_restart
  !----------------------------------------------------------------------------
  !
  ! ... this module contains methods to read and write data saved by the
  !     ballistic conductance code pwcond.x to restart smoothly
  !
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : tmp_dir, create_directory
  USE io_global, ONLY : ionode, ionode_id
  USE mp_images, ONLY : intra_image_comm
  USE mp,        ONLY : mp_bcast
  USE cond_files,ONLY : tran_prefix, tk_file
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: cond_writefile, cond_readfile
  !
  INTEGER, PRIVATE :: iunout
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE cond_writefile( what, kcurr, ecurr, tcurr )
      !------------------------------------------------------------------------
      !
      USE cond,                 ONLY : nenergy, earr, nkpts, xyk, wkpt
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER, INTENT(IN), OPTIONAL :: ecurr, kcurr
      REAL(DP), INTENT(IN), OPTIONAL :: tcurr
      !
      CHARACTER(LEN=256)  :: dirname, filename
      INTEGER :: ierr, ik
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      INTEGER, external :: find_free_unit

      ! look for an empty unit for transmission files,
      IF ( ionode )  iunout = find_free_unit( )
      !
      dirname = TRIM(tmp_dir) // TRIM(tran_prefix)
      !
      ! open the restart file
      IF ( ionode ) THEN
         !
         IF ( what=='init' ) THEN
            filename = TRIM(dirname) // '.rec'
         ELSEIF ( what=='tran' ) THEN
            filename = TRIM(dirname) // '_' // tk_file // '_k' // &
                 TRIM(int_to_char(kcurr)) // '_e' // TRIM(int_to_char(ecurr))
         ELSE
            CALL errore('cond_writefile','unknown what',1)
         ENDIF
         open (iunout, FILE=filename, FORM='formatted', IOSTAT=ierr)
         !
      END IF
      !
      CALL mp_bcast(ierr, ionode_id, intra_image_comm)
      !
      CALL errore('cond_writefile ', &
                  'cannot open recover file for writing', ierr )
      !
      IF ( ionode ) THEN
         !
         ! here we start writing the cond-punch-file
         IF ( what=='init' ) THEN
            !
            WRITE(iunout,"('NUMBER_OF_ENERGIES')")
            WRITE(iunout,*) nenergy
            WRITE(iunout,"('ENERGY_LIST')")
            WRITE(iunout,*) earr(:)
            !
            WRITE(iunout,"('NUMBER_OF_K-POINTS')")
            WRITE(iunout,*) nkpts
            DO ik = 1, nkpts
               !
               WRITE(iunout,"('K-POINT',i4)") ik
               WRITE(iunout,*) xyk(:,ik)
               WRITE(iunout,*) wkpt(ik)
               !
            END DO
            !
         ELSEIF ( what=='tran' ) THEN
            !
            WRITE(iunout,"('PARTIAL_TRANSMISSION')")
            WRITE(iunout,*) ecurr, kcurr, tcurr
            !
         ENDIF
         !
         CLOSE (unit=iunout, status='keep')
         !
      ENDIF
      RETURN
      !
    END SUBROUTINE cond_writefile
    !
    !------------------------------------------------------------------------
    SUBROUTINE cond_readfile( what, ierr, kcurr, ecurr, tcurr )
      !------------------------------------------------------------------------
      !
      USE cond,                 ONLY : nenergy, earr, nkpts, xyk, wkpt
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER, INTENT(IN), OPTIONAL :: ecurr, kcurr
      REAL(DP), INTENT(OUT), OPTIONAL :: tcurr
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=256)  :: dirname, filename
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      INTEGER, external :: find_free_unit
      INTEGER :: ne, ie, nk, ik
      LOGICAL :: exist
      REAL(DP) :: elist(nenergy), kpt(2), wk
      !
      ierr = 0
      !
      ! look for an empty unit for transmission files
      IF (ionode)  iunout = find_free_unit( )
      !
      dirname = TRIM(tmp_dir) // TRIM(tran_prefix)
      !
      IF ( what=='init' ) THEN
         filename = TRIM(dirname) // '.rec'
      ELSEIF ( what=='tran' ) THEN
         filename = TRIM(dirname) // '_' // tk_file // '_k' // &
              TRIM(int_to_char(kcurr)) // '_e' // TRIM(int_to_char(ecurr))
      ELSE
         CALL errore('cond_readfile','unknown what',1)
      ENDIF
      !
      INQUIRE (FILE=filename, EXIST=exist )
      IF ( .NOT.exist) CALL errore('cond_readfile','file not found',1)
      OPEN (iunout, FILE=filename, FORM='formatted', STATUS='old', IOSTAT=ierr)

      IF ( what=='init' ) THEN
         !
         READ(iunout,*)
         READ(iunout,*) ne
         IF ( ne .NE. nenergy ) ierr = 1
         READ(iunout,*)
         READ(iunout,*) elist(:)
         DO ie=1,ne
            IF (abs(elist(ie) - earr(ie)) .GT. 1.d-10) THEN
               ierr = ie+1
               EXIT
            ENDIF
         ENDDO
         !
         IF (ierr .NE. 0 ) CALL errore('cond_readfile', &
            'error while reading energies from info file',ierr)
         !
         READ(iunout,*)
         READ(iunout,*) nk
         IF ( nk .NE. nkpts ) ierr = 1
         !
         DO ik = 1, nk
            READ(iunout,*)
            READ(iunout,*) kpt(:)
            IF ( sum(abs(kpt(:) - xyk(:,ik))) .GT. 3.d-10 ) THEN
               ierr = ik+1
               EXIT
            ENDIF
            READ(iunout,*) wk
            IF ( abs(wk - wkpt(ik)) .GT. 1.d-10 ) THEN
               ierr = nk+ik+1
               EXIT
            ENDIF
         END DO
         IF (ierr .NE. 0 ) CALL errore('cond_readfile', &
              'error while reading k-points from info file',ierr)
         !
      ELSE
         !
         READ(iunout,*)
         READ(iunout,*) tcurr
         !
      END IF
      !
      CLOSE (iunout, STATUS='keep' )
      !
      RETURN
      !
    END SUBROUTINE cond_readfile
    !
END MODULE cond_restart
