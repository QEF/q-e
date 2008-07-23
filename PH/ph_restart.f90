!
! Copyright (C) 2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
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
  USE io_files,  ONLY : tmp_dir, prefix, iunpun, xmlpun, &
                        qexml_version, qexml_version_init
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
  PUBLIC :: ph_writefile, ph_readfile
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
  !
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE ph_writefile( what, irr )
      !------------------------------------------------------------------------
      !
      USE global_version,       ONLY : version_number
      USE control_ph,           ONLY : current_iq, xml_not_of_pw, done_bands, &
                                       ldisp, lnscf, epsil, trans, elph, zue
      USE ramanm,               ONLY : lraman, elop
      USE disp, ONLY : nqs, x_q, done_iq
      USE xml_io_base, ONLY : write_header, write_status_ph, &
                              write_control_ph, write_q, create_directory
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER, INTENT(IN) :: irr
      !
      CHARACTER(LEN=256)  :: dirname, filename
      INTEGER :: ierr
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
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
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.phsave'
      !
      ! ... create the main restart directory
      !
      CALL create_directory( dirname )
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
                    & TRIM( xmlpun ) // '.' // TRIM(int_to_char(current_iq))
            IF (what=='data') &
               CALL iotk_open_write( iunpun, FILE = TRIM( filename ), &
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
            CALL write_status_ph(current_iq, xml_not_of_pw, done_bands)
!-------------------------------------------------------------------------------
! ... CONTROL 
!-------------------------------------------------------------------------------
         !
            CALL write_control_ph( ldisp, lnscf, epsil, trans, elph, zue, &
                      lraman, elop )
         !
!-------------------------------------------------------------------------------
! ... Q AND K POINTS
!------------------------------------------------------------------------------
         !
            CALL write_q( nqs, x_q, done_iq )
         !
            CALL iotk_close_write( iunpun )
         ELSEIF (what=='data') THEN 
         !
!-------------------------------------------------------------------------------
! ... PARTIAL ITEMS
!-------------------------------------------------------------------------------
         !
             CALL write_partial_ph()
         !
         !
             CALL iotk_close_write( iunpun )
         ELSEIF (what=='data_dyn') THEN 

             CALL write_ph_dyn(filename,irr)

         END IF

      END IF


      RETURN
      !
      CONTAINS
      
         SUBROUTINE write_partial_ph()
            USE modes, ONLY : nirr, npert, u
            USE partial, ONLY : done_irr
            USE control_ph, ONLY : current_iq, epsil, trans, elph, zue, lgamma,&
                                   where_rec, rec_code
            USE ramanm,  ONLY : lraman, elop, ramtns, eloptns
            USE efield_mod, ONLY : zstareu, zstarue0, epsilon

            IMPLICIT NONE
            INTEGER :: imode0, imode, irr, ipert, iq, iunout
            CHARACTER(LEN=256) :: filename, filename1

            CALL iotk_write_begin( iunpun, "PARTIAL_PH" )
            !
            CALL iotk_write_dat(iunpun,"STOPPED_IN",where_rec)
            !
            CALL iotk_write_dat(iunpun,"RECOVER_CODE",rec_code)
            !
            CALL iotk_write_dat(iunpun,"QPOINT_NUMBER",current_iq)
            !
            IF (trans) THEN
              !
               CALL iotk_write_dat(iunpun,"NUMBER_IRR_REP",nirr)
              !
               imode0=0
               DO irr=1,nirr
                  CALL iotk_write_dat(iunpun,"NUMBER_OF_PERTURBATIONS",&
                                         npert(irr))
                  DO ipert=1,npert(irr)
                     imode=imode0+ipert
                     CALL iotk_write_dat(iunpun,"DISPLACEMENT_PATTERN",&
                                         u(:,imode))
                  ENDDO
                  imode0=imode0+npert(irr)
               ENDDO
            ENDIF
            IF (epsil.AND.lgamma) THEN
               CALL iotk_write_dat(iunpun,"DIELECTRIC_CONSTANT",epsilon)
               CALL iotk_write_dat(iunpun,"EFFECTIVE_CHARGES_EU",zstareu)
            ENDIF
            IF (zue.AND.lgamma) &
               CALL iotk_write_dat(iunpun,"EFFECTIVE_CHARGES_UE",zstarue0)
            IF (lraman.AND.lgamma) &
               CALL iotk_write_dat(iunpun,"RAMAN_TNS",ramtns)
            IF (elop.AND.lgamma) &
               CALL iotk_write_dat(iunpun,"ELOP_TNS",eloptns)
            !
             CALL iotk_write_end(iunpun, "PARTIAL_PH" )
            RETURN
         END SUBROUTINE write_partial_ph

        SUBROUTINE write_ph_dyn(filename, irr)
           USE partial, ONLY : done_irr
           USE dynmat,  ONLY : dyn_rec
           USE control_ph, ONLY : trans

           IMPLICIT NONE
           INTEGER :: irr, iunout
           CHARACTER(LEN=256) :: filename, filename1
           CHARACTER(LEN=6), EXTERNAL :: int_to_char

           IF (trans) THEN
              IF (done_irr(irr)) THEN
                 !
                 CALL iotk_free_unit( iunout, ierr )
                 !
                 filename1= TRIM(filename) // "." // TRIM(int_to_char(irr))
                 !
                 CALL iotk_open_write(iunout, FILE = TRIM(filename1), &
                                   BINARY = .FALSE., IERR = ierr )

                 CALL iotk_write_begin(iunout, "PARTIAL_MATRIX")
                 CALL iotk_write_dat(iunout, "DONE_IRR", done_irr(irr))
                 CALL iotk_write_dat(iunout, "PARTIAL_DYN", dyn_rec(:,:))
                 CALL iotk_write_end(iunout, "PARTIAL_MATRIX")
                 CALL iotk_close_write(iunout)
              ENDIF
           ENDIF
           RETURN
        END SUBROUTINE write_ph_dyn

    END SUBROUTINE ph_writefile 
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
      dirname = TRIM( tmp_dir ) // TRIM( prefix ) // '.phsave'
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
      USE control_ph, ONLY : current_iq, xml_not_of_pw, done_bands
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
         CALL iotk_scan_dat( iunpun, "XML_NOT_OF_PW", xml_not_of_pw )         
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
      CALL mp_bcast( xml_not_of_pw,    ionode_id, intra_image_comm )
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
      USE control_ph, ONLY : ldisp, lnscf, epsil, trans, elph, zue
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
         CALL iotk_scan_dat( iunpun, "BANDS_REQUIRED", lnscf )
         CALL iotk_scan_dat( iunpun, "ELECTRIC_FIELD", epsil )
         CALL iotk_scan_dat( iunpun, "PHONON_RUN", trans )
         CALL iotk_scan_dat( iunpun, "ELECTRON_PHONON", elph )
         CALL iotk_scan_dat( iunpun, "EFFECTIVE_CHARGE_PH", zue )
         CALL iotk_scan_dat( iunpun, "RAMAN_TENSOR", lraman )
         CALL iotk_scan_dat( iunpun, "ELECTRO_OPTIC", elop )
         !
         CALL iotk_scan_end( iunpun, "CONTROL" )
         !
         CALL iotk_close_read( iunpun )
      END IF
      CALL mp_bcast( ldisp,  ionode_id, intra_image_comm )
      CALL mp_bcast( lnscf,   ionode_id, intra_image_comm )
      CALL mp_bcast( epsil,   ionode_id, intra_image_comm )
      CALL mp_bcast( trans,   ionode_id, intra_image_comm )
      CALL mp_bcast( elph,    ionode_id, intra_image_comm )
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
         ALLOCATE(done_iq(nqs))
         CALL iotk_scan_dat( iunpun, "Q-POINT_COORDINATES", x_q(1:3,1:nqs) )
         !
         CALL iotk_scan_dat( iunpun, "Q-POINT_DONE", done_iq(1:nqs) )
         !
         CALL iotk_scan_end( iunpun, "Q_POINTS" )
         !
         CALL iotk_close_read( iunpun )
      ENDIF

      CALL mp_bcast( nqs,    ionode_id, intra_image_comm )

      IF (.NOT. ionode) THEN
         ALLOCATE(x_q(3,nqs))
         ALLOCATE(done_iq(nqs))
      ENDIF

      CALL mp_bcast( x_q,    ionode_id, intra_image_comm )
      CALL mp_bcast( done_iq,    ionode_id, intra_image_comm )

      RETURN
      !
    END SUBROUTINE read_q

    SUBROUTINE read_partial_ph( dirname, ierr )

    USE modes, ONLY : nirr
    USE partial, ONLY : done_irr, comp_irr
    USE disp, ONLY : done_iq
    USE dynmat,  ONLY : dyn_rec, dyn
    USE control_ph, ONLY : current_iq, trans

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: dirname
    INTEGER,          INTENT(OUT) :: ierr
    INTEGER :: irr, iunout

    CHARACTER(LEN=256)    :: filename, filename1
    CHARACTER(LEN=6), EXTERNAL :: int_to_char

    IF ( ionode ) THEN
       !
       filename= TRIM( dirname ) // '/' // &
               & TRIM( xmlpun ) // '.' // TRIM(int_to_char(current_iq))

       CALL iotk_open_read( iunpun, FILE = TRIM( filename ), IERR = ierr )
       !
    END IF
    !
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    !
    IF ( ierr > 0 ) RETURN
    !
    IF (ionode) THEN
       IF (trans) THEN
          done_irr=0
          dyn=(0.0_DP,0.0_DP)
          dyn_rec=(0.0_DP,0.0_DP)
          DO irr=0,nirr
             CALL iotk_free_unit( iunout, ierr )
             filename1=TRIM(filename) // "." // TRIM(int_to_char(irr))
             CALL iotk_open_read(iunout, FILE = TRIM(filename1), &
                                       BINARY = .FALSE., IERR = ierr )
             
             IF (ierr == 0 ) then
                CALL iotk_scan_begin( iunout, "PARTIAL_MATRIX" )
                CALL iotk_scan_dat(iunout,"DONE_IRR",done_irr(irr))
                IF (done_irr(irr)==1) comp_irr(irr)=1
                CALL iotk_scan_dat(iunout,"PARTIAL_DYN",&
                                                dyn_rec(:,:))
                dyn(:,:)=dyn(:,:) + dyn_rec(:,:)
                CALL iotk_scan_end( iunout, "PARTIAL_MATRIX" )
                CALL iotk_close_read( iunout )
             ELSE
                done_iq(current_iq)=0
                ierr=0
             END IF
          ENDDO
       ENDIF
    ENDIF
    IF (trans) THEN
       CALL mp_bcast( done_irr,  ionode_id, intra_image_comm )
       CALL mp_bcast( comp_irr,  ionode_id, intra_image_comm )
       CALL mp_bcast( dyn_rec,  ionode_id, intra_image_comm )
       CALL mp_bcast( dyn,  ionode_id, intra_image_comm )
    ENDIF

    RETURN
    END SUBROUTINE read_partial_ph

    SUBROUTINE read_u( dirname, ierr )

    USE modes, ONLY : nirr, npert, u
    USE control_ph, ONLY : current_iq, epsil, trans, elph, zue, lgamma, &
                           where_rec, rec_code
    USE ramanm,  ONLY : lraman, elop, ramtns, eloptns
    USE efield_mod, ONLY : zstareu, zstarue0, epsilon

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: dirname
    INTEGER,          INTENT(OUT) :: ierr
    INTEGER :: imode0, imode, irr, ipert, iq, iunout

    CHARACTER(LEN=256)    :: filename, filename1
    CHARACTER(LEN=6), EXTERNAL :: int_to_char

    IF ( ionode ) THEN
       !
       filename= TRIM( dirname ) // '/' // &
               & TRIM( xmlpun ) // '.' // TRIM(int_to_char(current_iq))

       CALL iotk_open_read( iunpun, FILE = TRIM( filename ), IERR = ierr )
       !
    END IF
    !
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    !
    IF ( ierr > 0 ) RETURN
    !
    IF (ionode) THEN
       CALL iotk_scan_begin( iunpun, "PARTIAL_PH" )
       !
       CALL iotk_scan_dat(iunpun,"STOPPED_IN",where_rec)
       !
       CALL iotk_scan_dat(iunpun,"RECOVER_CODE",rec_code)
       !
       CALL iotk_scan_dat(iunpun,"QPOINT_NUMBER",iq)
       IF (iq /= current_iq) CALL errore('read_partial_ph', &
              'problems with current_iq', 1 )
       IF (trans) THEN
          CALL iotk_scan_dat(iunpun,"NUMBER_IRR_REP",nirr)
          imode0=0
          DO irr=0,nirr
             IF (irr > 0) THEN
                CALL iotk_scan_dat(iunpun,"NUMBER_OF_PERTURBATIONS", npert(irr))
                DO ipert=1,npert(irr)
                   imode=imode0+ipert
                   CALL iotk_scan_dat(iunpun,"DISPLACEMENT_PATTERN",u(:,imode))
                ENDDO
                imode0=imode0+npert(irr)
             ENDIF
          ENDDO
       ENDIF
       IF (epsil.and.lgamma) THEN
          CALL iotk_scan_dat(iunpun,"DIELECTRIC_CONSTANT",epsilon)
          CALL iotk_scan_dat(iunpun,"EFFECTIVE_CHARGES_EU",zstareu)
       ENDIF
       IF (zue.and.lgamma) &
          CALL iotk_scan_dat(iunpun,"EFFECTIVE_CHARGES_UE",zstarue0)
       IF (lraman.and.lgamma) &
          CALL iotk_scan_dat(iunpun,"RAMAN_TNS",ramtns)
       IF (elop) &
          CALL iotk_scan_dat(iunpun,"ELOP_TNS",eloptns)
       CALL iotk_scan_end( iunpun, "PARTIAL_PH" )
       !
       CALL iotk_close_read( iunpun )
    ENDIF

    CALL mp_bcast( where_rec,  ionode_id, intra_image_comm )
    CALL mp_bcast( rec_code,  ionode_id, intra_image_comm )
    IF (trans) THEN
       CALL mp_bcast( nirr,  ionode_id, intra_image_comm )
       CALL mp_bcast( npert,  ionode_id, intra_image_comm )
       CALL mp_bcast( u,  ionode_id, intra_image_comm )
    ENDIF

    IF (epsil.and.lgamma) THEN
       CALL mp_bcast( epsilon, ionode_id, intra_image_comm )
       CALL mp_bcast( zstareu, ionode_id, intra_image_comm )
    ENDIF
    IF (zue.and.lgamma) CALL mp_bcast( zstarue0, ionode_id, intra_image_comm )
    IF (lraman.and.lgamma) CALL mp_bcast( ramtns, ionode_id, intra_image_comm )
    IF (elop.and.lgamma) CALL mp_bcast( eloptns,  ionode_id, intra_image_comm )

    RETURN
    END SUBROUTINE read_u
    !
END MODULE ph_restart
