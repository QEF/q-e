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
  USE iotk_module
  !
  USE kinds,     ONLY : DP
  USE qexml_module, ONLY : qexml_write_header
  USE xml_io_base, ONLY : create_directory, attr
  USE io_files,  ONLY : tmp_dir, xmlpun, iunpun, qexml_version, &
       qexml_version_init
  USE io_global, ONLY : ionode, ionode_id
  USE mp_global, ONLY : intra_image_comm
  USE mp,        ONLY : mp_bcast
  USE cond_files,  ONLY : tran_prefix, tk_file
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
  ! variables to describe qexml current version
  ! and back compatibility
  !
  LOGICAL :: qexml_version_before_1_4_0 = .FALSE.
  !
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE cond_writefile( what, kcurr, ecurr, tcurr )
      !------------------------------------------------------------------------
      !
      USE global_version,       ONLY : version_number
      USE cond,                 ONLY : nenergy, earr, nkpts, xyk, wkpt
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER, INTENT(IN), OPTIONAL :: ecurr, kcurr
      REAL(DP), INTENT(IN), OPTIONAL :: tcurr
      !
      CHARACTER(LEN=256)  :: dirname, filename
      INTEGER :: ierr
      CHARACTER(LEN=6), EXTERNAL :: int_to_char


      ! look for an empty unit for transmission files,
      ! (while info file goes in iunpun defined in io_files)
      IF ( ionode )  CALL iotk_free_unit(iunout, ierr)
      !
      CALL mp_bcast(ierr, ionode_id, intra_image_comm)
      !
      CALL errore('cond_writefile ', 'no free units to write ', ierr)
      !
      dirname = TRIM(tmp_dir) // TRIM(tran_prefix) // '.cond_save'
      !
      ! create the main restart directory
      CALL create_directory(dirname)
      !
      ! open the restart file
      IF ( ionode ) THEN
         !
         ! open XML descriptor
         ierr=0
         IF ( what=='init' ) THEN

            CALL iotk_open_write(iunpun, FILE=TRIM(dirname) // '/' // &
                                 TRIM(xmlpun), BINARY=.FALSE., IERR=ierr)
         ELSEIF ( what=='tran' ) THEN
            filename = TRIM(dirname) // '/' // tk_file // '_k' // &
               TRIM(int_to_char(kcurr)) // '_e' // TRIM(int_to_char(ecurr))
            CALL iotk_open_write(iunout, FILE=TRIM(filename), &
                                  BINARY=.FALSE., IERR=ierr)
         ELSE
            CALL errore('cond_writefile','unknown what',1)
         ENDIF
         !
      END IF
      !
      CALL mp_bcast(ierr, ionode_id, intra_image_comm)
      !
      CALL errore('cond_writefile ', &
                  'cannot open xml_recover file for writing', ierr )
      !
      IF ( ionode ) THEN
         !
         ! here we start writing the cond-punch-file
         IF ( what=='init' ) THEN
            !
            CALL qexml_write_header("PWCOND", TRIM(version_number))
            !
            CALL write_elist(nenergy, earr)
            !
            CALL write_klist(nkpts, xyk, wkpt)
            !
            CALL iotk_close_write(iunpun)
            !
         ELSEIF ( what=='tran' ) THEN
            !
            CALL write_transmission(ecurr, kcurr, tcurr)
            !
            CALL iotk_close_write(iunout)
            !
         ENDIF
         !
      ENDIF


      RETURN
      !
    END SUBROUTINE cond_writefile
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_elist( ne, elist )
      !------------------------------------------------------------------------
      !
      INTEGER,  INTENT(IN) :: ne
      REAL(DP), INTENT(IN) :: elist(:)
      !
      !
      CALL iotk_write_begin( iunpun, "SCATTERING_ENERGIES" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_ENERGIES", ne )
      !
      CALL iotk_write_attr( attr, "UNITS", "eV", FIRST = .TRUE. )
      !
      CALL iotk_write_dat( iunpun, "ENERGY_LIST", elist(:), ATTR=attr, COLUMNS=1 )
      !
      CALL iotk_write_end( iunpun, "SCATTERING_ENERGIES" )
      !
    END SUBROUTINE write_elist
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_klist( nk, klist, wk )
      !------------------------------------------------------------------------
      !
      INTEGER,  INTENT(IN) :: nk
      REAL(DP), INTENT(IN) :: klist(:,:), wk(:)
      !
      INTEGER :: ik
      !
      CALL iotk_write_begin( iunpun, "K-POINTS_MESH" )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      CALL iotk_write_empty( iunpun, "UNITS_FOR_K-POINTS", attr )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_K-POINTS", nk )
      !
      DO ik = 1, nk
         !
         CALL iotk_write_attr( attr, "XY", klist(:,ik), FIRST = .TRUE. )
         !
         CALL iotk_write_attr( attr, "WEIGHT", wk(ik) )
         !
         CALL iotk_write_empty( iunpun, "K-POINT" // &
                              & TRIM( iotk_index(ik) ), attr )
         !
      END DO
      !
      CALL iotk_write_end( iunpun, "K-POINTS_MESH" )
      !
    END SUBROUTINE write_klist
    !
    !------------------------------------------------------------------------
    SUBROUTINE write_transmission( ie, ik, t )
      !------------------------------------------------------------------------
      !
      INTEGER,  INTENT(IN) :: ie, ik
      REAL(DP), INTENT(IN) :: t
      !
      CALL iotk_write_dat( iunout, "PARTIAL_TRANSMISSION", t )
      !
    END SUBROUTINE write_transmission
    !
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
      !
      CHARACTER(LEN=256) :: dirname
      !
      ierr = 0
      !
      dirname = TRIM( tmp_dir ) // TRIM( tran_prefix ) // '.cond_save'
      !
      ! look for an empty unit for transmission files
      IF (ionode)  CALL iotk_free_unit( iunout, ierr )
      !
      CALL mp_bcast( ierr,ionode_id,intra_image_comm )
      !
      CALL errore( 'cond_readfile', &
                   'no free units to read restart file', ierr )
      !
      SELECT CASE( what )
      CASE( 'init' )
         !
         qexml_version_init = .FALSE.
         CALL read_header( dirname, ierr )
         IF (ierr .NE. 0 ) CALL errore('cond_readfile', &
            'error while reading header of info file',ierr)
         !
         CALL read_elist( dirname, nenergy, earr, ierr )
         IF (ierr .NE. 0 ) CALL errore('cond_readfile', &
            'error while reading energies from info file',ierr)
         !
         CALL read_klist( dirname, nkpts, xyk, wkpt, ierr )
         IF (ierr .NE. 0 ) CALL errore('cond_readfile', &
            'error while reading k-points from info file',ierr)
         !
      CASE( 'tran' )
         !
         CALL read_transmission( dirname, kcurr, ecurr, tcurr, ierr )
         ! the corresponding file may not be present for all (e,k)
         !
      END SELECT
      !
      RETURN
      !
    END SUBROUTINE cond_readfile
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_header( dirname, ierr )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the format version of the current xml datafile
      !
      USE parser, ONLY : version_compare
      USE xml_io_base
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr

      ierr = 0
      IF ( qexml_version_init ) RETURN
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr.NE.0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "HEADER" )
         !
         CALL iotk_scan_empty( iunpun, "FORMAT", ATTR=attr )
         !
         CALL iotk_scan_attr( attr, "VERSION", qexml_version )
         !
         qexml_version_init = .TRUE.
         !
         CALL iotk_scan_end( iunpun, "HEADER" )
         !
         CALL iotk_close_read( iunpun )
         !
      ENDIF
      !
      CALL mp_bcast( qexml_version,       ionode_id, intra_image_comm )
      CALL mp_bcast( qexml_version_init,  ionode_id, intra_image_comm )

      !
      ! init logical variables for versioning
      !
      qexml_version_before_1_4_0 = .FALSE.
      !
      IF ( TRIM( version_compare( qexml_version, "1.4.0" )) == "older" ) &
         qexml_version_before_1_4_0 = .TRUE.
      !
      RETURN
    END SUBROUTINE read_header
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_elist( dirname, ne, elist, ierr )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,  INTENT(IN) :: ne
      REAL(DP), INTENT(IN) :: elist(:)
      INTEGER,  INTENT(OUT) :: ierr
      ! local
      INTEGER :: ne_, ie
      REAL(DP) :: elist_(ne)
      !
      ierr = 0
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr.NE.0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "SCATTERING_ENERGIES" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_ENERGIES", ne_ )
         IF ( ne_.NE.ne ) ierr = 1
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      IF ( ierr .NE. 0 ) RETURN
      !
      IF ( ionode ) THEN
         CALL iotk_scan_dat( iunpun, "ENERGY_LIST", elist_ )
         DO ie=1,ne
            IF (abs(elist_(ie) - elist(ie)) .GT. 1.d-10) THEN
               ierr = ie+1
               EXIT
            ENDIF
         ENDDO
         !
         CALL iotk_scan_end( iunpun, "SCATTERING_ENERGIES" )
         !
         CALL iotk_close_read( iunpun )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
    END SUBROUTINE read_elist
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_klist( dirname, nk, klist, wk, ierr )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,  INTENT(IN) :: nk
      REAL(DP), INTENT(IN) :: klist(:,:), wk(:)
      INTEGER,  INTENT(OUT) :: ierr
      !
      INTEGER :: nk_, ik
      REAL(DP) :: kpt_(2), wk_
      !
      ierr = 0
      !
      IF ( ionode ) &
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr .NE. 0 ) RETURN
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "K-POINTS_MESH" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_K-POINTS", nk_ )
         !
         IF ( nk_ .NE. nk ) ierr = 1
         !
      ENDIF

      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr .NE. 0 ) RETURN

      IF ( ionode ) THEN
         DO ik = 1, nk
            !
            CALL iotk_scan_empty(iunpun, "K-POINT"//TRIM(iotk_index(ik)), attr)
            !
            CALL iotk_scan_attr( attr, "XY", kpt_ )
            IF ( sum(abs(kpt_(:) - klist(:,ik))) .GT. 3.d-10 ) THEN
               ierr = ik+1
               EXIT
            ENDIF
            !
            CALL iotk_scan_attr( attr, "WEIGHT", wk_ )
            !
            IF ( abs(wk_ - wk(ik)) .GT. 1.d-10 ) THEN
               ierr = nk+ik+1
               EXIT
            ENDIF
            !
         END DO
         !
         CALL iotk_scan_end( iunpun, "K-POINTS_MESH" )
         !
         CALL iotk_close_read( iunpun )
         !
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
    END SUBROUTINE read_klist
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_transmission( dirname, ik, ie, t, ierr )
      !------------------------------------------------------------------------
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER, INTENT(IN) :: ie, ik
      REAL(DP), INTENT(OUT) :: t
      INTEGER,  INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=256) :: filename
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      ierr = 0
      !
      IF ( ionode ) THEN
      !
         filename = TRIM( dirname ) // '/' // tk_file // '_k' // &
            TRIM(int_to_char(ik)) // '_e' // TRIM(int_to_char(ie))
         CALL iotk_open_read( iunout, FILE = TRIM( filename ), &
                               BINARY = .FALSE., IERR = ierr )
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr .NE. 0 ) THEN
          ierr = 1 ! file not found
          RETURN
      ENDIF
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_dat( iunout, "PARTIAL_TRANSMISSION", t, IERR = ierr )
         !
         CALL iotk_close_read( iunout )
         !
      ENDIF
      !
      IF ( ierr .NE. 0 )  ierr = 2 ! file not readable?
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
    END SUBROUTINE read_transmission
    !
    !
    !
END MODULE cond_restart
