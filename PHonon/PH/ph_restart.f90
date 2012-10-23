!
! Copyright (C) 2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE ph_restart
  !----------------------------------------------------------------------------
  !
  ! ... this module contains methods to read and write data saved by the
  !     phonon code to restart smoothly
  !
  USE iotk_module
  !
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : prefix, iunpun, xmlpun, xmlpun_base, &
                        qexml_version, qexml_version_init
  USE control_ph, ONLY : tmp_dir_ph
  USE io_global, ONLY : ionode, ionode_id
  USE mp_global, ONLY : intra_image_comm
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: ph_writefile, ph_readfile, init_status_run, check_status_run, &
            destroy_status_run
  !
  INTEGER, PRIVATE :: iunout
  !
  LOGICAL :: lheader           = .FALSE., &
             lstatus_read      = .FALSE., &
             lcontrol_ph_read  = .FALSE., &
             lq_read           = .FALSE., &
             lu_read           = .FALSE., &
             lpartial_dyn_read   = .FALSE.
  !
  ! variables to describe qexml current version
  ! and back compatibility
  !
  LOGICAL :: qexml_version_before_1_4_0 = .FALSE.

  CHARACTER(iotk_attlenx)  :: attr
  !
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE ph_writefile( what, irr )
      !------------------------------------------------------------------------
      !
      USE global_version,       ONLY : version_number
      USE control_ph,           ONLY : current_iq, done_bands, &
                                       ldisp, epsil, trans, zue, zeu
      USE ramanm,               ONLY : lraman, elop
      USE el_phon,              ONLY : elph
      USE disp, ONLY : nqs, x_q
      USE xml_io_base, ONLY : write_header, create_directory
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER, INTENT(IN) :: irr
      !
      CHARACTER(LEN=256)  :: dirname, filename
      INTEGER :: ierr
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      LOGICAL :: exst
      !
      IF ( ionode ) THEN
         !
         ! ... look for an empty unit (only ionode needs it)
         !
         CALL iotk_free_unit( iunpun, ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'ph_writefile ', &
                   'no free units to write ', ierr )
      !
      dirname = TRIM( tmp_dir_ph ) // TRIM( prefix ) // '.phsave'
      !
      ! ... create the main restart directory if it is not there
      !
      IF (ionode) inquire (file =TRIM(dirname)//'/data-file.xml', &
                           exist = exst)
      !
      CALL mp_bcast( exst, ionode_id, intra_image_comm )
      !
      IF (.NOT. exst) CALL create_directory( dirname )
      !
      ! ... open the ph_recover file
      !
      IF ( ionode ) THEN
         !
         ! ... open XML descriptor
         !
         ierr=0
         IF (what=='init') THEN
            CALL iotk_open_write( iunpun, FILE = TRIM( dirname ) // '/' // &
                             & TRIM( xmlpun ), BINARY = .FALSE., IERR = ierr )
         ELSEIF (what=='data' .OR. what=='data_dyn') THEN
            filename= TRIM( dirname ) // '/' // &
                    & TRIM( xmlpun_base ) // '.' // TRIM(int_to_char(current_iq))
            IF (what=='data') &
               CALL iotk_open_write( iunpun, FILE = TRIM( filename ) // '.xml', &
                                          BINARY = .FALSE., IERR = ierr )
         ELSE
            CALL errore('ph_writefile','unknown what',1)
         ENDIF
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'ph_writefile ', &
                   'cannot open xml_recover file for writing', ierr )
      !
      IF ( ionode ) THEN
         !
         ! ... here we start writing the ph-punch-file
         !
!-------------------------------------------------------------------------------
! ... HEADER
!-------------------------------------------------------------------------------
         !
         IF (what=='init') THEN
         !
            CALL write_header( "PH", TRIM(version_number) )
         !
!-------------------------------------------------------------------------------
! ... STATUS
!-------------------------------------------------------------------------------
!   With what='init' the routine writes the status of the code
!   needed to control the main flow of the dispersion calculation.
!   The mesh of q points, the current q, the main flags that control
!   the run and the flag done_bands, that tells if this routine is
!   called before or after the band calculation.
!
            CALL write_status_ph(current_iq, done_bands)
!-------------------------------------------------------------------------------
! ... CONTROL
!-------------------------------------------------------------------------------
         !
            CALL write_control_ph( ldisp, epsil, trans, elph, zue, zeu, &
                      lraman, elop )
         !
!-------------------------------------------------------------------------------
! ... Q POINTS
!------------------------------------------------------------------------------
         !
            CALL write_q( nqs, x_q )
         !
            CALL iotk_close_write( iunpun )
         ELSEIF (what=='data') THEN
         !
!-------------------------------------------------------------------------------
! ... PARTIAL ITEMS
!-------------------------------------------------------------------------------
!
! with what='data' this routine writes the information on the irreducible
! representations. Number of irreducible representations, number of modes
! for each representation and displacements u. Moreover it writes the
! point in which the phonon code has written the recover file: where_rec
! and rec_code give the same information, where_rec is a string
! and rec_code an integer. The former is easy to read in the xml file,
! the latter is simpler to use in the code. Moreover this routine saves
! the tensors that contain the result of the calculations done so far:
! epsilon, zstareu, ramtns, eloptns, dyn, zstarue.
!
             CALL write_partial_ph()
         !
         !
             CALL iotk_close_write( iunpun )
         ELSEIF (what=='data_dyn') THEN
!
! with what='data_dyn' this routine writes the information calculated
! separately for each irreducible representation. The contributions
! of the representation to the dynamical matrix and to the Born effective
! charges dP/du.
!
             CALL write_ph_dyn(filename,irr)

         END IF

      END IF


      RETURN
      !
      CONTAINS

         SUBROUTINE write_partial_ph()
            USE modes, ONLY : nirr, npert, u, name_rap_mode, nsymq
            USE partial, ONLY : done_irr
            USE control_ph, ONLY : current_iq, epsil, trans, zue, lgamma,&
                                   where_rec, rec_code, done_epsil, done_zeu, &
                                   done_zue
            USE el_phon, ONLY : elph
            USE ramanm,  ONLY : lraman, elop, ramtns, eloptns, done_lraman, &
                                done_elop
            USE efield_mod, ONLY : zstareu, zstarue, epsilon

            IMPLICIT NONE
            INTEGER :: imode0, imode, irr, ipert, iq

            CALL iotk_write_begin( iunpun, "TENSOR_INFO" )
            !
            CALL iotk_write_dat(iunpun,"STOPPED_IN",where_rec)
            !
            CALL iotk_write_dat(iunpun,"RECOVER_CODE",rec_code)
            !
            CALL iotk_write_dat(iunpun,"QPOINT_NUMBER",current_iq)
            !
            CALL iotk_write_dat(iunpun,"QPOINT_GROUP_RANK",nsymq)
            !
!            IF (trans.OR.zeu.OR.elph) THEN
              !
               CALL iotk_write_dat(iunpun,"NUMBER_IRR_REP",nirr)
              !
               imode0=0
               DO irr=1,nirr
                  CALL iotk_write_dat(iunpun,"NUMBER_OF_PERTURBATIONS",&
                                         npert(irr))
                  CALL iotk_write_dat(iunpun,"SYMMETRY_TYPE",&
                                         name_rap_mode(imode0+1))
                  DO ipert=1,npert(irr)
                     imode=imode0+ipert
                     CALL iotk_write_dat(iunpun,"DISPLACEMENT_PATTERN",&
                                         u(:,imode))
                  ENDDO
                  imode0=imode0+npert(irr)
               ENDDO
!            ENDIF
!
!   Save the current flags
!
            IF (lgamma) THEN
               CALL iotk_write_dat( iunpun,"DONE_ELECTRIC_FIELD",done_epsil )
               CALL iotk_write_dat( iunpun,"DONE_EFFECTIVE_CHARGE_EU",done_zeu )
               CALL iotk_write_dat( iunpun,"DONE_EFFECTIVE_CHARGE_PH",done_zue )
               CALL iotk_write_dat( iunpun,"DONE_RAMAN_TENSOR",done_lraman )
               CALL iotk_write_dat( iunpun,"DONE_ELECTRO_OPTIC",done_elop )
!
!    save all calculated tensors
!
               IF (done_epsil) &
                  CALL iotk_write_dat(iunpun,"DIELECTRIC_CONSTANT",epsilon)
               IF (done_zeu) &
                  CALL iotk_write_dat(iunpun,"EFFECTIVE_CHARGES_EU",zstareu)
               IF (done_lraman) &
                  CALL iotk_write_dat(iunpun,"RAMAN_TNS",ramtns)
               IF (done_elop) &
                  CALL iotk_write_dat(iunpun,"ELOP_TNS",eloptns)
               IF (done_zue) &
                  CALL iotk_write_dat(iunpun,"EFFECTIVE_CHARGES_UE",zstarue)
            ENDIF
            !
             CALL iotk_write_end(iunpun, "TENSOR_INFO" )
            RETURN
         END SUBROUTINE write_partial_ph

        SUBROUTINE write_ph_dyn(filename, irr)
           USE partial, ONLY : done_irr
           USE dynmat,  ONLY : dyn_rec
           USE efield_mod, ONLY : zstarue0_rec
           USE control_ph, ONLY : trans, zue

           IMPLICIT NONE
           INTEGER :: irr, iunout
           CHARACTER(LEN=256) :: filename, filename1
           CHARACTER(LEN=6), EXTERNAL :: int_to_char

           IF (trans) THEN
              IF (done_irr(irr)/=0) THEN
                 !
                 CALL iotk_free_unit( iunout, ierr )
                 !
                 filename1= TRIM(filename) // "." // TRIM(int_to_char(irr))
                 !
                 CALL iotk_open_write(iunout, FILE = TRIM(filename1)//'.xml', &
                                   BINARY = .FALSE., IERR = ierr )

                 CALL iotk_write_begin(iunout, "PARTIAL_MATRIX")
                 CALL iotk_write_dat(iunout, "DONE_IRR", done_irr(irr))
                 CALL iotk_write_dat(iunout, "PARTIAL_DYN", dyn_rec(:,:))
                 IF (zue) CALL iotk_write_dat(iunout, &
                                           "PARTIAL_ZUE", zstarue0_rec(:,:))
                 CALL iotk_write_end(iunout, "PARTIAL_MATRIX")
                 CALL iotk_close_write(iunout)
              ENDIF
           ENDIF
           RETURN
        END SUBROUTINE write_ph_dyn

    END SUBROUTINE ph_writefile

    SUBROUTINE write_control_ph( ldisp, epsil, trans, elph, zue, zeu, &
                      lraman, elop)
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: ldisp, epsil, trans, elph, zue, zeu, &
                      lraman, elop

      CALL iotk_write_begin( iunpun, "CONTROL" )
      !
      CALL iotk_write_dat( iunpun, "DISPERSION_RUN", ldisp )
      CALL iotk_write_dat( iunpun, "ELECTRIC_FIELD", epsil )
      CALL iotk_write_dat( iunpun, "PHONON_RUN", trans )
      CALL iotk_write_dat( iunpun, "ELECTRON_PHONON", elph )
      CALL iotk_write_dat( iunpun, "EFFECTIVE_CHARGE_EU", zeu )
      CALL iotk_write_dat( iunpun, "EFFECTIVE_CHARGE_PH", zue )
      CALL iotk_write_dat( iunpun, "RAMAN_TENSOR", lraman )
      CALL iotk_write_dat( iunpun, "ELECTRO_OPTIC", elop )
      !
      CALL iotk_write_end( iunpun, "CONTROL" )
      !
      RETURN
    END SUBROUTINE write_control_ph

    SUBROUTINE write_status_ph(current_iq, done_bands)
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: current_iq
      LOGICAL, INTENT(IN) :: done_bands

      CALL iotk_write_begin( iunpun, "STATUS_PH" )
      !
      CALL iotk_write_dat( iunpun, "DONE_BANDS", done_bands )
      CALL iotk_write_dat( iunpun, "CURRENT_Q", current_iq )
      !
      CALL iotk_write_end( iunpun, "STATUS_PH" )
      !
      RETURN
    END SUBROUTINE write_status_ph
    !

    SUBROUTINE write_q( nqs, x_q)
      !------------------------------------------------------------------------
      !
      INTEGER, INTENT(IN) :: nqs
      REAL(DP), INTENT(IN) :: x_q(3,nqs)
      !
      CALL iotk_write_begin( iunpun, "Q_POINTS" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_Q_POINTS", nqs  )
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      !
      CALL iotk_write_empty( iunpun, "UNITS_FOR_Q-POINT", attr )
      !
      CALL iotk_write_dat( iunpun, "Q-POINT_COORDINATES", x_q(:,:), COLUMNS=3 )
      !
      CALL iotk_write_end( iunpun, "Q_POINTS" )
      !
      RETURN
    END SUBROUTINE write_q
    !
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE ph_readfile( what, ierr )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: what
      INTEGER,          INTENT(OUT) :: ierr
      !
      CHARACTER(LEN=256) :: dirname
      !
      ierr = 0
      !
      dirname = TRIM( tmp_dir_ph ) // TRIM( prefix ) // '.phsave'
      !
      ! ... look for an empty unit
      !
      IF (ionode) THEN
         CALL iotk_free_unit( iunout, ierr )
      ENDIF
      CALL mp_bcast(ierr,ionode_id,intra_image_comm)
      !
      CALL errore( 'ph_readfile', &
                   'no free units to read wavefunctions', ierr )
      !
      lheader = .FALSE.
      lstatus_read      = .FALSE.
      lcontrol_ph_read  = .FALSE.
      lq_read          = .FALSE.
      lpartial_dyn_read   = .FALSE.
      lu_read   = .FALSE.
      !
      SELECT CASE( what )
      CASE( 'init' )
         !
         lheader = .TRUE.
         lstatus_read=.TRUE.
         lq_read=.TRUE.
         lcontrol_ph_read=.TRUE.
         !
      CASE( 'data' )
         !
         lpartial_dyn_read = .TRUE.
         !
      CASE( 'data_u' )
         !
         lu_read = .TRUE.
         !
      CASE( 'reset' )
         !
         lheader = .FALSE.
         lstatus_read      = .FALSE.
         lcontrol_ph_read  = .FALSE.
         lq_read           = .FALSE.
         lpartial_dyn_read   = .FALSE.
         lu_read   = .FALSE.
         !
         RETURN
         !
      END SELECT
      !
      IF ( lheader ) THEN
         !
         CALL read_header( dirname, ierr )
         IF (ierr /= 0 ) RETURN
         !
      ENDIF
      IF ( lstatus_read ) THEN
         !
         CALL read_status_ph( dirname, ierr )
         IF ( ierr /= 0 ) RETURN
         !
      END IF
      IF ( lcontrol_ph_read ) THEN
         !
         CALL read_control_ph( dirname, ierr )
         IF ( ierr /= 0 ) RETURN
         !
      END IF
      IF ( lq_read ) THEN
         !
         CALL read_q( dirname, ierr )
         IF ( ierr /= 0 ) RETURN
         !
      END IF
      IF ( lpartial_dyn_read ) THEN
         !
         CALL read_partial_ph( dirname, ierr )
         IF ( ierr /= 0 ) RETURN
         !
      END IF
      IF ( lu_read ) THEN
         !
         CALL read_u( dirname, ierr )
         IF ( ierr /= 0 ) RETURN
         !
      END IF
      !
      RETURN
      !
    END SUBROUTINE ph_readfile
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_header( dirname, ierr )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the format version of the current xml datafile
      !
      USE parser, ONLY : version_compare
      USE xml_io_base, ONLY : attr

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
      IF ( ierr /=0 ) RETURN
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
    !------------------------------------------------------------------------
    SUBROUTINE read_status_ph( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE control_ph, ONLY : current_iq, done_bands
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      !
      IF ( ionode ) THEN
         !
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      ! ... then selected tags are read from the other sections
      !
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "STATUS_PH" )
         !
         CALL iotk_scan_dat( iunpun, "DONE_BANDS", done_bands )
         !
         CALL iotk_scan_dat( iunpun, "CURRENT_Q", current_iq )
         !
         CALL iotk_scan_end( iunpun, "STATUS_PH" )
         !
         CALL iotk_close_read( iunpun )
      END IF
      !
      CALL mp_bcast( done_bands,       ionode_id, intra_image_comm )
      CALL mp_bcast( current_iq,       ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_status_ph
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_control_ph( dirname, ierr )
      !------------------------------------------------------------------------
      USE control_ph, ONLY : ldisp, epsil, trans, zue, zeu
      USE el_phon, ONLY : elph
      USE ramanm,     ONLY : lraman, elop

      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      IF ( ionode ) THEN
         !
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF ( ionode ) THEN
         CALL iotk_scan_begin( iunpun, "CONTROL" )
         !
         CALL iotk_scan_dat( iunpun, "DISPERSION_RUN", ldisp )
         CALL iotk_scan_dat( iunpun, "ELECTRIC_FIELD", epsil )
         CALL iotk_scan_dat( iunpun, "PHONON_RUN", trans )
         CALL iotk_scan_dat( iunpun, "ELECTRON_PHONON", elph )
         CALL iotk_scan_dat( iunpun, "EFFECTIVE_CHARGE_EU", zeu )
         CALL iotk_scan_dat( iunpun, "EFFECTIVE_CHARGE_PH", zue )
         CALL iotk_scan_dat( iunpun, "RAMAN_TENSOR", lraman )
         CALL iotk_scan_dat( iunpun, "ELECTRO_OPTIC", elop )
         !
         CALL iotk_scan_end( iunpun, "CONTROL" )
         !
         CALL iotk_close_read( iunpun )
      END IF
      CALL mp_bcast( ldisp,  ionode_id, intra_image_comm )
      CALL mp_bcast( epsil,   ionode_id, intra_image_comm )
      CALL mp_bcast( trans,   ionode_id, intra_image_comm )
      CALL mp_bcast( elph,    ionode_id, intra_image_comm )
      CALL mp_bcast( zeu,     ionode_id, intra_image_comm )
      CALL mp_bcast( zue,     ionode_id, intra_image_comm )
      CALL mp_bcast( lraman,  ionode_id, intra_image_comm )
      CALL mp_bcast( elop,    ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_control_ph
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_q( dirname, ierr )
      !------------------------------------------------------------------------
      !
      USE disp, ONLY : nqs, x_q, done_iq
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: dirname
      INTEGER,          INTENT(OUT) :: ierr
      !
      IF ( ionode ) THEN
         !
         CALL iotk_open_read( iunpun, FILE = TRIM( dirname ) // '/' // &
                            & TRIM( xmlpun ), IERR = ierr )
         !
      END IF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      IF ( ierr > 0 ) RETURN
      !
      IF (ionode) THEN
         CALL iotk_scan_begin( iunpun, "Q_POINTS" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_Q_POINTS", nqs  )
         !
         ALLOCATE(x_q(3,nqs))
         CALL iotk_scan_dat( iunpun, "Q-POINT_COORDINATES", x_q(1:3,1:nqs) )
         !
         CALL iotk_scan_end( iunpun, "Q_POINTS" )
         !
         CALL iotk_close_read( iunpun )
      ENDIF

      CALL mp_bcast( nqs,    ionode_id, intra_image_comm )

      IF (.NOT. ionode) THEN
         ALLOCATE(x_q(3,nqs))
      ENDIF

      CALL mp_bcast( x_q,    ionode_id, intra_image_comm )

      RETURN
      !
    END SUBROUTINE read_q

    SUBROUTINE read_partial_ph( dirname, ierr )

    USE modes, ONLY : nirr
    USE partial, ONLY : done_irr, comp_irr
    USE efield_mod, ONLY : zstarue0_rec, zstarue0
    USE dynmat,  ONLY : dyn_rec, dyn
    USE control_ph, ONLY : current_iq, trans, zue

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: dirname
    INTEGER,          INTENT(OUT) :: ierr
    INTEGER :: irr, iunout

    CHARACTER(LEN=256)    :: filename, filename1
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    IF ( ionode ) THEN
       !
       filename= TRIM( dirname ) // '/' // &
               & TRIM( xmlpun_base ) // '.' // TRIM(int_to_char(current_iq))
       !
    END IF

    IF (ionode) THEN
       IF (trans) THEN
          done_irr=0
          dyn=(0.0_DP,0.0_DP)
          dyn_rec=(0.0_DP,0.0_DP)
          zstarue0=(0.0_DP, 0.0_DP)
          DO irr=0,nirr
             CALL iotk_free_unit( iunout, ierr )
             filename1=TRIM(filename) // "." // TRIM(int_to_char(irr))
             CALL iotk_open_read(iunout, FILE = TRIM(filename1) // '.xml', &
                                       BINARY = .FALSE., IERR = ierr )

             IF (ierr == 0 ) then
                CALL iotk_scan_begin( iunout, "PARTIAL_MATRIX" )
                CALL iotk_scan_dat(iunout,"DONE_IRR",done_irr(irr))
                IF (done_irr(irr)==1) comp_irr(irr)=1
                CALL iotk_scan_dat(iunout,"PARTIAL_DYN",&
                                                dyn_rec(:,:))
                dyn(:,:)=dyn(:,:) + dyn_rec(:,:)
                IF (zue) THEN
                   CALL iotk_scan_dat(iunout, &
                                     "PARTIAL_ZUE", zstarue0_rec(:,:))
                   zstarue0(:,:)=zstarue0(:,:)+zstarue0_rec(:,:)
                ENDIF
                CALL iotk_scan_end( iunout, "PARTIAL_MATRIX" )
                CALL iotk_close_read( iunout )
             ELSE
                ierr=0
             END IF
          ENDDO
       ENDIF
    ENDIF
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    IF (trans) THEN
       CALL mp_bcast( done_irr,  ionode_id, intra_image_comm )
       CALL mp_bcast( comp_irr,  ionode_id, intra_image_comm )
       CALL mp_bcast( dyn_rec,  ionode_id, intra_image_comm )
       CALL mp_bcast( dyn,  ionode_id, intra_image_comm )
       IF (zue) THEN
          CALL mp_bcast( zstarue0,  ionode_id, intra_image_comm )
          CALL mp_bcast( zstarue0_rec,  ionode_id, intra_image_comm )
       ENDIF
    ENDIF

    RETURN
    END SUBROUTINE read_partial_ph

    SUBROUTINE read_u( dirname, ierr )
!
!   This routine reads the tensors that have been already calculated and
!   the displacement patterns. It reads also the code that tells the phonon
!   where it stopped. The convention is the following:
!
!   rec_code    where_rec     status description
!
!    -1000                    Nothing has been read. There is no recover file.
!    -40         phq_setup    Only the displacements u have been read from file
!    -30         phq_init     u and dyn(0) read from file
!    -25                      not yet active. Restart in solve_e_fpol
!    -20         solve_e      all previous. Stopped within solve_e. There
!                             should be a recover file.
!    -10         solve_e2     epsilon and zstareu are available if requested.
!                             Within solve_e2. There should be a recover file.
!     2          phescf       all previous, raman tenson and elop tensor are
!                             available if required.
!     10         solve_linter all previous, within solve linter. Recover file
!                             should be present.
!     20         phqscf       all previous dyn_rec(irr) and zstarue0(irr) are
!                             available.
!     30         dynmatrix    all previous, dyn and zstarue are available.
!
!
    USE modes, ONLY : nirr, npert, u, name_rap_mode
    USE control_ph, ONLY : current_iq, epsil, trans, zue, zeu, lgamma, &
                           where_rec, rec_code, rec_code_read, done_epsil, &
                           done_zeu, done_zue
    USE el_phon, ONLY : elph
    USE ramanm,  ONLY : lraman, elop, ramtns, eloptns, done_lraman, done_elop
    USE efield_mod, ONLY : zstareu, zstarue, epsilon

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: dirname
    INTEGER,          INTENT(OUT) :: ierr
    INTEGER :: imode0, imode, irr, ipert, iq, iunout

    CHARACTER(LEN=256)    :: filename
    CHARACTER(LEN=6), EXTERNAL :: int_to_char

    IF ( ionode ) THEN
       !
       filename= TRIM( dirname ) // '/' // &
               & TRIM( xmlpun_base ) // '.' // TRIM(int_to_char(current_iq))

       CALL iotk_open_read( iunpun, FILE = TRIM( filename ) // '.xml', IERR = ierr )
       !
    END IF
    !
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    !
    IF ( ierr > 0 ) RETURN
    !
    IF (ionode) THEN
       CALL iotk_scan_begin( iunpun, "TENSOR_INFO" )
       !
       CALL iotk_scan_dat(iunpun,"STOPPED_IN",where_rec)
       !
       CALL iotk_scan_dat(iunpun,"RECOVER_CODE",rec_code_read)
       !
       CALL iotk_scan_dat(iunpun,"QPOINT_NUMBER",iq)
       !
    ENDIF
    CALL mp_bcast( iq,  ionode_id, intra_image_comm )
    IF (iq /= current_iq) CALL errore('read_partial_ph', &
              'problems with current_iq', 1 )

    IF (ionode) THEN
!       IF (trans.OR.zeu.OR.elph) THEN
          CALL iotk_scan_dat(iunpun,"NUMBER_IRR_REP",nirr)
          imode0=0
          DO irr=0,nirr
             IF (irr > 0) THEN
                CALL iotk_scan_dat(iunpun,"NUMBER_OF_PERTURBATIONS", npert(irr))
                CALL iotk_scan_dat(iunpun,"SYMMETRY_TYPE",&
                                                    name_rap_mode(imode0+1))
                DO ipert=1,npert(irr)
                   imode=imode0+ipert
                   CALL iotk_scan_dat(iunpun,"DISPLACEMENT_PATTERN",u(:,imode))
                   IF (ipert > 1) name_rap_mode(imode)=name_rap_mode(imode0+1)
                ENDDO
                imode0=imode0+npert(irr)
             ENDIF
          ENDDO
!       ENDIF
!
!   read all flags
!
       IF (lgamma) THEN
          CALL iotk_scan_dat( iunpun, "DONE_ELECTRIC_FIELD", done_epsil )
          CALL iotk_scan_dat( iunpun, "DONE_EFFECTIVE_CHARGE_EU", done_zeu )
          CALL iotk_scan_dat( iunpun, "DONE_EFFECTIVE_CHARGE_PH", done_zue )
          CALL iotk_scan_dat( iunpun, "DONE_RAMAN_TENSOR", done_lraman )
          CALL iotk_scan_dat( iunpun, "DONE_ELECTRO_OPTIC", done_elop )

          IF (done_epsil) &
             CALL iotk_scan_dat(iunpun,"DIELECTRIC_CONSTANT",epsilon)
          IF (done_zeu) &
             CALL iotk_scan_dat(iunpun,"EFFECTIVE_CHARGES_EU",zstareu)
          IF (done_lraman) &
             CALL iotk_scan_dat(iunpun,"RAMAN_TNS",ramtns)
          IF (done_elop) &
             CALL iotk_scan_dat(iunpun,"ELOP_TNS",eloptns)
          IF (done_zue) &
             CALL iotk_scan_dat(iunpun,"EFFECTIVE_CHARGES_UE",zstarue)
       ENDIF
       !
       CALL iotk_scan_end( iunpun, "TENSOR_INFO" )
       !
       CALL iotk_close_read( iunpun )
    ENDIF

    CALL mp_bcast( where_rec,  ionode_id, intra_image_comm )
    CALL mp_bcast( rec_code_read,  ionode_id, intra_image_comm )
    rec_code=rec_code_read
    IF (lgamma) THEN
       CALL mp_bcast( done_epsil,  ionode_id, intra_image_comm )
       CALL mp_bcast( done_zeu,  ionode_id, intra_image_comm )
       CALL mp_bcast( done_zue,  ionode_id, intra_image_comm )
       CALL mp_bcast( done_lraman,  ionode_id, intra_image_comm )
       CALL mp_bcast( done_elop,  ionode_id, intra_image_comm )
       IF (done_epsil) CALL mp_bcast( epsilon, ionode_id, intra_image_comm )
       IF (done_zeu) CALL mp_bcast( zstareu, ionode_id, intra_image_comm )
       IF (done_zue) CALL mp_bcast( zstarue, ionode_id, intra_image_comm )
       IF (done_lraman) CALL mp_bcast( ramtns, ionode_id, intra_image_comm )
       IF (done_elop) CALL mp_bcast( eloptns,  ionode_id, intra_image_comm )
    ENDIF

!    IF (trans.OR.zeu.OR.elph) THEN
       CALL mp_bcast( nirr,  ionode_id, intra_image_comm )
       CALL mp_bcast( npert,  ionode_id, intra_image_comm )
       CALL mp_bcast( u,  ionode_id, intra_image_comm )
       CALL mp_bcast( name_rap_mode,  ionode_id, intra_image_comm )
!    ENDIF

    RETURN
    END SUBROUTINE read_u

  !----------------------------------------------------------------------------
    SUBROUTINE check_status_run(  )
  !----------------------------------------------------------------------------
  !
  ! ...
  ! ... This routine sets the situation of the grid according to
  ! ... the files that it finds on the directory .phsave.
  ! ... Check if representation files exist and which representations
  ! ... have been already calculated.
  ! ... set the initial information on the grid
  ! ... it sets done_iq and done_rep_iq to 1 for the q and the
  ! ... representations that have already been done.
  ! ... Moreover it sets rep_iq, the number of representation for each q.
  !
  USE disp, ONLY : nqs, done_iq, done_rep_iq, rep_iq, nsymq_iq, npert_iq
  USE ions_base, ONLY : nat
  USE control_ph, ONLY : trans, zeu
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256)  :: dirname, filename, filename1
  INTEGER :: iunout, iq, irr, ierr
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  CALL init_status_run()
  !
  dirname = TRIM( tmp_dir_ph ) // TRIM( prefix ) // '.phsave'
  DO iq=1, nqs
     IF ( ionode ) THEN
        !
        CALL iotk_free_unit( iunout, ierr )
        filename= TRIM( dirname ) // '/' // &
                & TRIM( xmlpun_base ) // '.' // TRIM(int_to_char(iq))

        CALL iotk_open_read( iunout, FILE = TRIM( filename ) // '.xml', IERR = ierr )

        IF (ierr /= 0) CYCLE
        CALL iotk_scan_begin( iunout, "TENSOR_INFO" )
        !
        CALL iotk_scan_dat(iunout,"QPOINT_GROUP_RANK",nsymq_iq(iq))
        !
        IF (trans.OR.zeu) THEN
           CALL iotk_scan_dat(iunout,"NUMBER_IRR_REP",rep_iq(iq))
           DO irr=1,rep_iq(iq)
              CALL iotk_scan_dat(iunout,"NUMBER_OF_PERTURBATIONS",&
                                         npert_iq(irr,iq))
           ENDDO
        ENDIF

        CALL iotk_scan_end( iunout, "TENSOR_INFO" )
        CALL iotk_close_read(iunout)

        IF (trans.OR.zeu) THEN
           DO irr=rep_iq(iq)+1,3*nat
              done_rep_iq(irr,iq)=2
           ENDDO
           DO irr=0,rep_iq(iq)
              filename1=TRIM(filename) // "." // TRIM(int_to_char(irr))
              CALL iotk_open_read(iunout, FILE = TRIM(filename1) // '.xml', &
                                       BINARY = .FALSE., IERR = ierr )
              IF (ierr /= 0 ) CYCLE
              CALL iotk_scan_begin( iunout, "PARTIAL_MATRIX" )
              CALL iotk_scan_dat(iunout,"DONE_IRR",done_rep_iq(irr,iq))
              CALL iotk_scan_end( iunout, "PARTIAL_MATRIX" )
              CALL iotk_close_read(iunout)
           END DO
           done_iq(iq)=1
           DO irr=0,rep_iq(iq)
              IF (done_rep_iq(irr,iq) /= 1) THEN
                 done_iq(iq)=0
                 EXIT
              ENDIF
           ENDDO
        END IF
     END IF
  END DO
  !
  CALL mp_bcast( done_iq, ionode_id, intra_image_comm )
  CALL mp_bcast( rep_iq, ionode_id, intra_image_comm )
  CALL mp_bcast( nsymq_iq, ionode_id, intra_image_comm )
  CALL mp_bcast( npert_iq, ionode_id, intra_image_comm )
  CALL mp_bcast( done_rep_iq, ionode_id, intra_image_comm )
  !
  !
  RETURN
  !
  END SUBROUTINE check_status_run


   SUBROUTINE init_status_run()
   USE  disp, ONLY : nqs, done_iq, done_rep_iq, rep_iq, nsymq_iq, npert_iq, &
                     comp_iq, comp_irr_iq
   USE ions_base, ONLY : nat
   USE mp_global, ONLY : nimage
   IMPLICIT NONE

   ALLOCATE(done_iq(nqs))
   ALLOCATE(rep_iq(nqs))
   ALLOCATE(done_rep_iq(0:3*nat,nqs))
   ALLOCATE(nsymq_iq(nqs))
   ALLOCATE(npert_iq(3*nat,nqs))
   ALLOCATE(comp_iq(nqs))
   IF (nimage>1) ALLOCATE(comp_irr_iq(0:3*nat,nqs))

   done_iq=0
   rep_iq=3*nat
   done_rep_iq=0
   nsymq_iq=0
   npert_iq=0
   IF ( allocated(comp_irr_iq) ) comp_irr_iq=0
   comp_iq=0

   RETURN
   END SUBROUTINE init_status_run

   SUBROUTINE destroy_status_run()
   USE start_k,    ONLY : xk_start, wk_start
   USE disp, ONLY : nqs, done_iq, done_rep_iq, rep_iq, nsymq_iq, npert_iq, &
                    comp_iq, comp_irr_iq, x_q
   IMPLICIT NONE

   IF (ALLOCATED(x_q)) DEALLOCATE(x_q)
   IF (ALLOCATED(done_iq)) DEALLOCATE(done_iq)
   IF (ALLOCATED(rep_iq)) DEALLOCATE(rep_iq)
   IF (ALLOCATED(done_rep_iq)) DEALLOCATE(done_rep_iq)
   IF (ALLOCATED(nsymq_iq)) DEALLOCATE(nsymq_iq)
   IF (ALLOCATED(npert_iq)) DEALLOCATE(npert_iq)
   IF (ALLOCATED(comp_iq)) DEALLOCATE(comp_iq)
   IF (ALLOCATED(comp_irr_iq)) DEALLOCATE(comp_irr_iq)
!
! Note that these two variables are allocated by read_file.
! They cannot be deallocated by clean_pw because the starting xk and wk
! points must be known at each q point.
! The logic of these two variables must be improved.
!
  IF (ALLOCATED( xk_start )) DEALLOCATE( xk_start )
  IF (ALLOCATED( wk_start )) DEALLOCATE( wk_start )

   END SUBROUTINE destroy_status_run

    !
END MODULE ph_restart
