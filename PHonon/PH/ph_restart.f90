!
! Copyright (C) 2008-2012 Quantum ESPRESSO group
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
  USE io_files,  ONLY : prefix, xmlpun, qexml_version, qexml_version_init
  USE control_ph, ONLY : tmp_dir_ph
  USE io_global, ONLY : ionode, ionode_id
  USE mp_images, ONLY : intra_image_comm
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: ph_writefile, ph_readfile, allocate_grid_variables, &
                         check_directory_phsave, destroy_status_run, &
                         check_available_bands
  !
  INTEGER :: iunpun
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
    SUBROUTINE ph_writefile( what, iq, irr, ierr )
      !------------------------------------------------------------------------
      !
      USE global_version,       ONLY : version_number
      USE control_ph,           ONLY : ldisp, epsil, trans, zue, zeu
      USE el_phon,              ONLY : elph
      USE freq_ph,              ONLY : fpol, nfs, fiu, current_iu
      USE ramanm,               ONLY : lraman, elop
      USE disp,                 ONLY : nqs, x_q, nq1, nq2, nq3
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: what
      INTEGER, INTENT(IN) :: iq, irr
      !
      INTEGER, INTENT(OUT) :: ierr
      !
      CALL ph_restart_set_filename( what, irr, iq, 1, ierr)
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
            CALL write_header_ph( "PH", TRIM(version_number) )
         !
!
!   With what='init' the routine writes the main variables that
!   control the main flow of the dispersion calculation:
!   The main flags of the phonon code, the mesh of q point, the
!   number of q points and their coordinates.
!
!-------------------------------------------------------------------------------
! ... CONTROL 
!-------------------------------------------------------------------------------
         !
            CALL write_control_ph( ldisp, epsil, trans, elph, zue, zeu, &
                      lraman, elop, fpol ) 
         !
!-------------------------------------------------------------------------------
! ... Q POINTS AND FREQUENCY POINTS
!------------------------------------------------------------------------------
         !
            CALL write_qu( nqs, nq1, nq2, nq3, x_q, nfs, fiu, fpol )
         !
         ELSEIF (what=='status_ph') THEN
!
!   In this case we save the information on the status of the calculation.
!   The current q point, the current frequency, the label and the
!   code with the point where the code arrived so far. 
!   The former is easy to read in the xml file, 
!   the latter is simpler to use in the code. 
!
            CALL write_status_ph(iq, current_iu)
            !
         ELSEIF (what=='data_u') THEN 
!
! with what='data_u' this routine writes the information on the irreducible
! representations. Number of irreducible representations, number of modes
! for each representation and displacements u.  
!
            CALL write_modes(iq)
         !
         ELSEIF (what=='polarization') THEN 
!
! With what='polarization' this routine saves the tensors that contain the 
! polarization as a function of frequency.
!
             CALL write_polarization(irr)
         !
         ELSEIF (what=='tensors') THEN 
!
! With what='tensors' this routine saves the tensors that contain the 
! result of the calculations done so far: epsilon, zstareu, ramtns, eloptns, 
! dyn, zstarue.
!
             CALL write_tensors()
         !
         ELSEIF (what=='data_dyn') THEN 
!
! with what='data_dyn' this routine writes the information calculated
! separately for each irreducible representation. The contributions
! of the representation to the dynamical matrix and to the Born effective 
! charges dP/du.
!
             CALL write_ph_dyn(irr)

         ELSEIF (what=='el_phon') THEN 
! with what='data_dyn' this routine writes the information calculated
! for this irreducible representation to the electron phonon    
!
             CALL write_el_phon(irr)

         END IF
         CALL iotk_close_write( iunpun )
      END IF


      RETURN
      !
      CONTAINS
      
         SUBROUTINE write_polarization(iu)
!
            USE freq_ph, ONLY : polar, done_iu, fpol, done_fpol, fiu

            IMPLICIT NONE
            INTEGER :: iu

            IF (.NOT.fpol) RETURN
            CALL iotk_write_begin( iunpun, "POLARIZ_IU" )
!
!   Save the current flags
!
            CALL iotk_write_dat( iunpun,"DONE_POLARIZ_IU",done_fpol )
!
!    Here we save the frequency dependent polarization at this iu
!
            CALL iotk_write_dat( iunpun, "FREQUENCY_IN_RY", fiu(iu) )
            CALL iotk_write_dat( iunpun, "CALCULATED_FREQUENCY", &
                                              done_iu(iu))
            IF ( done_iu(iu) ) &
                   CALL iotk_write_dat( iunpun, "POLARIZATION_IU",     &
                             polar(:,:,iu), COLUMNS=3 )
            !
            CALL iotk_write_end(iunpun, "POLARIZ_IU" )
            RETURN
         END SUBROUTINE write_polarization

         SUBROUTINE write_tensors()
!
            USE control_ph, ONLY : done_epsil, done_start_zstar, done_zeu, done_zue
            USE ramanm,  ONLY : lraman, elop, ramtns, eloptns, done_lraman, &
                                done_elop

            USE efield_mod, ONLY : zstareu0, zstareu, zstarue, epsilon

            IMPLICIT NONE

            CALL iotk_write_begin( iunpun, "EF_TENSORS" )
!
!   Save the current flags
!
            CALL iotk_write_dat( iunpun,"DONE_ELECTRIC_FIELD",done_epsil )
            CALL iotk_write_dat( iunpun,"DONE_START_EFFECTIVE_CHARGE",done_start_zstar )
            CALL iotk_write_dat( iunpun,"DONE_EFFECTIVE_CHARGE_EU",done_zeu )
            CALL iotk_write_dat( iunpun,"DONE_EFFECTIVE_CHARGE_PH",done_zue )
            CALL iotk_write_dat( iunpun,"DONE_RAMAN_TENSOR",done_lraman )
            CALL iotk_write_dat( iunpun,"DONE_ELECTRO_OPTIC",done_elop )
!
!    save all calculated tensors
!
            IF (done_epsil) &
               CALL iotk_write_dat(iunpun,"DIELECTRIC_CONSTANT",epsilon,COLUMNS=3)
            IF (done_start_zstar) &
               CALL iotk_write_dat(iunpun,"START_EFFECTIVE_CHARGES",zstareu0,COLUMNS=3)
            IF (done_zeu) &
               CALL iotk_write_dat(iunpun,"EFFECTIVE_CHARGES_EU",zstareu,COLUMNS=3)
            IF (done_lraman) &
               CALL iotk_write_dat(iunpun,"RAMAN_TNS",ramtns,COLUMNS=3)

            IF (done_elop) &
               CALL iotk_write_dat(iunpun,"ELOP_TNS",eloptns,COLUMNS=3)

            IF (done_zue) &
               CALL iotk_write_dat(iunpun,"EFFECTIVE_CHARGES_UE",zstarue)
            !
            CALL iotk_write_end(iunpun, "EF_TENSORS" )
            RETURN
         END SUBROUTINE write_tensors

         SUBROUTINE write_modes(iq)
            USE modes, ONLY : nirr, npert, u, name_rap_mode, num_rap_mode

            USE lr_symm_base, ONLY : nsymq, minus_q

            IMPLICIT NONE
            INTEGER :: imode0, imode, irr, ipert, iq

            CALL iotk_write_begin( iunpun, "IRREPS_INFO" )
            !
            CALL iotk_write_dat(iunpun,"QPOINT_NUMBER",iq)
            !
            CALL iotk_write_dat(iunpun,"QPOINT_GROUP_RANK",nsymq)
            !
            CALL iotk_write_dat(iunpun,"MINUS_Q_SYM",minus_q)
            !
            CALL iotk_write_dat(iunpun,"NUMBER_IRR_REP",nirr)
            !
            imode0=0
            DO irr=1,nirr
               CALL iotk_write_begin( iunpun, "REPRESENTION"// &
                                     TRIM( iotk_index( irr ) ) )
               CALL iotk_write_dat(iunpun,"NUMBER_OF_PERTURBATIONS",&
                                      npert(irr))
               DO ipert=1,npert(irr)
                  imode=imode0+ipert
                  CALL iotk_write_begin( iunpun, "PERTURBATION"// &
                                     TRIM( iotk_index( ipert ) ) )
                  CALL iotk_write_dat(iunpun,"SYMMETRY_TYPE_CODE", &
                                                      num_rap_mode(imode))

                  CALL iotk_write_dat(iunpun,"SYMMETRY_TYPE",&
                                       name_rap_mode(imode))
                  CALL iotk_write_dat(iunpun,"DISPLACEMENT_PATTERN",&
                                       u(:,imode))
                  CALL iotk_write_end( iunpun, "PERTURBATION"// &
                                     TRIM( iotk_index( ipert ) ) )
               ENDDO
               imode0=imode0+npert(irr)
               CALL iotk_write_end( iunpun, "REPRESENTION"// &
                                     TRIM( iotk_index( irr ) ) )
            ENDDO
            !
            CALL iotk_write_end(iunpun, "IRREPS_INFO" )
            RETURN
         END SUBROUTINE write_modes

        SUBROUTINE write_ph_dyn(irr)
           USE partial, ONLY : done_irr
           USE dynmat,  ONLY : dyn_rec
           USE efield_mod, ONLY : zstarue0_rec
           USE control_ph, ONLY : trans, zue

           IMPLICIT NONE
           INTEGER, INTENT(IN) :: irr

           IF (trans.OR.zeu) THEN
              IF (done_irr(irr)) THEN
                 !
                 CALL iotk_write_begin(iunpun, "PM_HEADER")
                 CALL iotk_write_dat(iunpun, "DONE_IRR", done_irr(irr))
                 CALL iotk_write_end(iunpun, "PM_HEADER")
                 CALL iotk_write_begin(iunpun, "PARTIAL_MATRIX")
                 CALL iotk_write_dat(iunpun, "PARTIAL_DYN", dyn_rec(:,:))
                 IF (zue.and.irr>0) CALL iotk_write_dat(iunpun, &
                                           "PARTIAL_ZUE", zstarue0_rec(:,:))

                 CALL iotk_write_end(iunpun, "PARTIAL_MATRIX")
              ENDIF
           ENDIF
           RETURN
        END SUBROUTINE write_ph_dyn

        SUBROUTINE write_el_phon(irr)
           USE el_phon, ONLY : done_elph, el_ph_mat_rec_col, elph
           USE klist, ONLY : nks
           USE wvfct, ONLY: nbnd
           USE qpoint, ONLY : nksqtot, xk_col
           USE control_lr, ONLY : lgamma
           IMPLICIT NONE
           INTEGER, INTENT(IN) :: irr
           INTEGER :: ik, ikk

           IF (.NOT. elph .OR. .NOT. done_elph(irr)) RETURN
           !
           CALL iotk_write_begin(iunpun, "EL_PHON_HEADER")
           CALL iotk_write_dat(iunpun, "DONE_ELPH", done_elph(irr))
           CALL iotk_write_end(iunpun, "EL_PHON_HEADER")
           CALL iotk_write_begin(iunpun, "PARTIAL_EL_PHON")
           CALL iotk_write_dat(iunpun, "NUMBER_OF_K", nksqtot)
           CALL iotk_write_dat(iunpun, "NUMBER_OF_BANDS", nbnd)
           DO ik=1,nksqtot
              ikk = 2 * ik - 1
              IF (lgamma) ikk = ik 
              CALL iotk_write_begin(iunpun, "K_POINT" // &
                                        TRIM( iotk_index( ik ) ) )
              CALL iotk_write_dat(iunpun, "COORDINATES_XK", &
                                         xk_col(:,ikk), COLUMNS=3)
              CALL iotk_write_dat(iunpun, "PARTIAL_ELPH", &
                                         el_ph_mat_rec_col(:,:,ik,:))
              CALL iotk_write_end(iunpun, "K_POINT" // &
                                        TRIM( iotk_index( ik ) ) )
           ENDDO
           CALL iotk_write_end(iunpun, "PARTIAL_EL_PHON")
        RETURN
        END SUBROUTINE write_el_phon

    END SUBROUTINE ph_writefile 

    !------------------------------------------------------------------------
    SUBROUTINE write_header_ph( creator_name, creator_version ) 
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(IN) :: creator_name, creator_version
      CHARACTER(5),  PARAMETER :: fmt_name = "QEXML"
      CHARACTER(5),  PARAMETER :: fmt_version = "1.4.0"

      CALL iotk_write_begin( iunpun, "HEADER" )
      !
      CALL iotk_write_attr(attr, "NAME",TRIM(fmt_name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "VERSION",TRIM(fmt_version) )
      CALL iotk_write_empty( iunpun, "FORMAT", ATTR=attr )
      !
      CALL iotk_write_attr(attr, "NAME",TRIM(creator_name), FIRST=.TRUE.)
      CALL iotk_write_attr(attr, "VERSION",TRIM(creator_version) )
      CALL iotk_write_empty( iunpun, "CREATOR", ATTR=attr )
      !
      CALL iotk_write_end( iunpun, "HEADER" )
      !
    END SUBROUTINE write_header_ph
    !
    !
    SUBROUTINE write_control_ph( ldisp, epsil, trans, elph, zue, zeu, &
                      lraman, elop, fpol) 
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: ldisp, epsil, trans, elph, zue, zeu, &
                      lraman, elop, fpol

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
      CALL iotk_write_dat( iunpun, "FREQUENCY_DEP_POL", fpol )
      !
      CALL iotk_write_end( iunpun, "CONTROL" )
      !
      RETURN
    END SUBROUTINE write_control_ph

    SUBROUTINE write_status_ph(current_iq, current_iu)
      !------------------------------------------------------------------------
      !
      USE control_ph, ONLY : where_rec, rec_code
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: current_iq, current_iu

      CALL iotk_write_begin( iunpun, "STATUS_PH" )
      CALL iotk_write_dat( iunpun, "STOPPED_IN", where_rec )
      CALL iotk_write_dat( iunpun, "RECOVER_CODE", rec_code )
      CALL iotk_write_dat( iunpun, "CURRENT_Q", current_iq )
      CALL iotk_write_dat( iunpun, "CURRENT_IU", current_iu )
      CALL iotk_write_end( iunpun, "STATUS_PH" )
      !
      RETURN
    END SUBROUTINE write_status_ph
    !

    SUBROUTINE write_qu( nqs, nq1, nq2, nq3, x_q, nfs, fiu, fpol)
      !------------------------------------------------------------------------
      !
      INTEGER, INTENT(IN) :: nqs, nfs, nq1, nq2, nq3
      REAL(DP), INTENT(IN) :: x_q(3,nqs), fiu(nfs)
      LOGICAL, INTENT(IN) :: fpol
      INTEGER :: nqind(3)
      !
      CALL iotk_write_begin( iunpun, "Q_POINTS" )
      !
      CALL iotk_write_dat( iunpun, "NUMBER_OF_Q_POINTS", nqs  )
      !
      IF (nqs > 1) THEN
         nqind(1) = nq1
         nqind(2) = nq2
         nqind(3) = nq3
         !
         CALL iotk_write_dat( iunpun, "MESH_DIMENSIONS", nqind(:), COLUMNS=3)
      ENDIF
      !
      CALL iotk_write_attr( attr, "UNITS", "2 pi / a", FIRST = .TRUE. )
      !
      CALL iotk_write_empty( iunpun, "UNITS_FOR_Q-POINT", attr )
      !
      CALL iotk_write_dat( iunpun, "Q-POINT_COORDINATES", x_q(:,:), COLUMNS=3 )
      !
      CALL iotk_write_end( iunpun, "Q_POINTS" )
      !
      IF (fpol) THEN
         !
         CALL iotk_write_begin( iunpun, "FREQUENCIES" )
         !
         CALL iotk_write_dat( iunpun, "NUMBER_OF_FREQUENCIES", nfs  )
         !
         CALL iotk_write_dat( iunpun, "FREQUENCY_VALUES", fiu(:), COLUMNS=1 )
         !
         CALL iotk_write_end( iunpun, "FREQUENCIES" )
         !
      ENDIF
      !
      RETURN
    END SUBROUTINE write_qu
    !
    !
    !------------------------------------------------------------------------
    SUBROUTINE ph_readfile( what, iq, irr, ierr )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN)  :: what
      INTEGER,          INTENT(IN)  :: irr, iq
      !
      !  irreducible representation and q point
      !
      INTEGER,          INTENT(OUT) :: ierr
      !
      CALL ph_restart_set_filename( what, irr, iq, -1, ierr)
      IF (ierr /= 0) RETURN
      !
      SELECT CASE( what )
      CASE( 'init' )
         !
         CALL read_header( ierr )
         IF (ierr /= 0 ) RETURN
         CALL read_control_ph( ierr )
         IF ( ierr /= 0 ) RETURN
         CALL read_qu( ierr )
         IF ( ierr /= 0 ) RETURN
         !
      CASE( 'status_ph')
         !
         CALL read_status_ph( ierr )
         IF ( ierr /= 0 ) RETURN
         !
      CASE( 'data_u' )
         !
         CALL read_modes( iq, ierr )
         IF ( ierr /= 0 ) RETURN
         !
      CASE( 'polarization' )
         !
         CALL read_polarization( irr, ierr )
         IF ( ierr /= 0 ) RETURN
         !
      CASE( 'tensors' )
         !
         CALL read_tensors( ierr )
         IF ( ierr /= 0 ) RETURN
         !
      CASE( 'data_dyn' )
         !
         CALL read_partial_ph( irr, ierr )
         IF ( ierr /= 0 ) RETURN
         !
      CASE( 'el_phon' )
         !
         CALL read_el_phon( irr, ierr )
         IF ( ierr /= 0 ) RETURN
         !
      CASE DEFAULT
         !
         CALL errore('ph_readfile','called with the wrong what',1)
         !
      END SELECT
      !
      IF (ionode) CALL iotk_close_read( iunpun )
      !
      RETURN
      !
    END SUBROUTINE ph_readfile
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_header( ierr )
      !------------------------------------------------------------------------
      !
      ! ... this routine reads the format version of the current xml datafile
      !
      USE parser, ONLY : version_compare
      USE xml_io_base, ONLY : attr

      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr

      ierr = 0
      IF ( qexml_version_init ) RETURN
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
    SUBROUTINE read_status_ph( ierr )
      !------------------------------------------------------------------------
      !
      !  This routine reads the status of ph. It tells where the code stopped
      !  There is both a number, to be used within the code, and a label
      !  that is easier to read within the recover file.
      !
      !   The convention is the following:
      !
      !  rec_code   where_rec     status description
      !
      !    -1000              Nothing has been read. There is no recover file.
      !    -40     phq_setup  Only the displacements u have been read from file
      !    -30     phq_init   u and dyn(0) read from file
      !    -25                not yet active. Restart in solve_e_fpol
      !    -20     solve_e    all previous. Stopped within solve_e. There 
      !                       should be a recover file.
      !    -10     solve_e2   epsilon and zstareu are available if requested. 
      !                       Within solve_e2. There should be a recover file.
      !     2      phescf     all previous, raman tenson and elop tensor are
      !                       available if required.
      !     10     solve_linter all previous, within solve linter. Recover file
      !                       should be present.
      !     20     phqscf     all previous dyn_rec(irr) and zstarue0(irr) are
      !                       available.
      !     30     dynmatrix  all previous, dyn and zstarue are available.
      !
      !
      USE control_ph, ONLY : current_iq, where_rec, rec_code_read
      USE freq_ph, ONLY : current_iu
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      !
      ! ... then selected tags are read from the other sections
      !
      ierr=0
      IF ( ionode ) THEN
         !
         CALL iotk_scan_begin( iunpun, "STATUS_PH" )
         CALL iotk_scan_dat( iunpun, "STOPPED_IN", where_rec )
         CALL iotk_scan_dat( iunpun, "RECOVER_CODE", rec_code_read )
         CALL iotk_scan_dat( iunpun, "CURRENT_Q", current_iq )         
         CALL iotk_scan_dat( iunpun, "CURRENT_IU", current_iu )         
         CALL iotk_scan_end( iunpun, "STATUS_PH" )
         !
      END IF
      !
      CALL mp_bcast( where_rec,        ionode_id, intra_image_comm )
      CALL mp_bcast( rec_code_read,         ionode_id, intra_image_comm )
      CALL mp_bcast( current_iq,       ionode_id, intra_image_comm )
      CALL mp_bcast( current_iu,       ionode_id, intra_image_comm )
      !
      RETURN
      !
    END SUBROUTINE read_status_ph
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_control_ph( ierr )
      !------------------------------------------------------------------------
      USE control_ph, ONLY : ldisp, epsil, trans, zue, zeu
      USE el_phon,    ONLY : elph
      USE ramanm,     ONLY : lraman, elop
      USE freq_ph,    ONLY : fpol
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      LOGICAL :: ldisp_, epsil_, trans_, zue_, zeu_, elph_, lraman_, elop_, &
                 fpol_
      !
      ierr=0
      IF ( ionode ) THEN
         CALL iotk_scan_begin( iunpun, "CONTROL" )
         !
         CALL iotk_scan_dat( iunpun, "DISPERSION_RUN", ldisp_ )
         CALL iotk_scan_dat( iunpun, "ELECTRIC_FIELD", epsil_ )
         CALL iotk_scan_dat( iunpun, "PHONON_RUN", trans_ )
         CALL iotk_scan_dat( iunpun, "ELECTRON_PHONON", elph_ )
         CALL iotk_scan_dat( iunpun, "EFFECTIVE_CHARGE_EU", zeu_ )
         CALL iotk_scan_dat( iunpun, "EFFECTIVE_CHARGE_PH", zue_ )
         CALL iotk_scan_dat( iunpun, "RAMAN_TENSOR", lraman_ )
         CALL iotk_scan_dat( iunpun, "ELECTRO_OPTIC", elop_ )
         CALL iotk_scan_dat( iunpun, "FREQUENCY_DEP_POL", fpol_ )
         !
         CALL iotk_scan_end( iunpun, "CONTROL" )
         !
      END IF
      CALL mp_bcast( ldisp_,   ionode_id, intra_image_comm )
      CALL mp_bcast( epsil_,   ionode_id, intra_image_comm )
      CALL mp_bcast( trans_,   ionode_id, intra_image_comm )
      CALL mp_bcast( elph_,    ionode_id, intra_image_comm )
      CALL mp_bcast( zeu_,     ionode_id, intra_image_comm )
      CALL mp_bcast( zue_,     ionode_id, intra_image_comm )
      CALL mp_bcast( lraman_,  ionode_id, intra_image_comm )
      CALL mp_bcast( elop_,    ionode_id, intra_image_comm )
      CALL mp_bcast( fpol_,    ionode_id, intra_image_comm )
      !
      IF (ldisp_ .neqv. ldisp) CALL errore('read_control_ph','wrong ldisp',1)
      IF (epsil_ .neqv. epsil) CALL errore('read_control_ph','wrong epsil',1)
      IF (trans_ .neqv. trans) CALL errore('read_control_ph','wrong trans',1)
      IF (elph_ .neqv. elph) CALL errore('read_control_ph','wrong elph',1)
      IF (zeu_ .neqv. zeu) CALL errore('read_control_ph','wrong zeu',1)
      IF (zue_ .neqv. zue) CALL errore('read_control_ph','wrong zue',1)
      IF (lraman_ .neqv. lraman) CALL errore('read_control_ph','wrong lraman',1)
      IF (elop_ .neqv. elop) CALL errore('read_control_ph','wrong elop',1)
      IF (fpol_ .neqv. fpol) CALL errore('read_control_ph','wrong fpol',1)
      !
      RETURN
      !
    END SUBROUTINE read_control_ph
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_qu( ierr )
      !------------------------------------------------------------------------
      !
      USE disp, ONLY : nqs, x_q, nq1, nq2, nq3, lgamma_iq
      USE freq_ph, ONLY : fpol, nfs, fiu
      !
      IMPLICIT NONE
      !
      INTEGER,          INTENT(OUT) :: ierr
      INTEGER :: nfs_, nqind(3), iq
      LOGICAL :: exst
      !
      ierr=0
      IF (ionode) THEN
         CALL iotk_scan_begin( iunpun, "Q_POINTS" )
         !
         CALL iotk_scan_dat( iunpun, "NUMBER_OF_Q_POINTS", nqs  )
         !
         IF (nqs > 1) & 
            CALL iotk_scan_dat( iunpun, "MESH_DIMENSIONS", nqind(1:3) )
         !
         ALLOCATE(x_q(3,nqs))
         CALL iotk_scan_dat( iunpun, "Q-POINT_COORDINATES", x_q(1:3,1:nqs) )
         !
         CALL iotk_scan_end( iunpun, "Q_POINTS" )
         !
         IF (fpol) THEN
            !
            CALL iotk_scan_begin( iunpun, "FREQUENCIES" )
            !
            CALL iotk_scan_dat( iunpun, "NUMBER_OF_FREQUENCIES", nfs_  )
            !
            CALL iotk_scan_dat( iunpun, "FREQUENCY_VALUES", fiu(1:nfs_) )
            !
            CALL iotk_scan_end( iunpun, "FREQUENCIES" )
            !
         ENDIF
      ENDIF

      CALL mp_bcast( nqs,    ionode_id, intra_image_comm )

      IF (nqs > 1) THEN
         CALL mp_bcast( nqind,    ionode_id, intra_image_comm )
         IF ( (nqind(1) /= nq1 ) .OR. (nqind(2) /= nq2) .OR. &
                 (nqind(3) /= nq3 ) )  &
               CALL errore('read_qu','nq1, nq2, or nq3 do not match',1)
            !
      ENDIF

      IF (.NOT. ionode) ALLOCATE(x_q(3,nqs))
      CALL mp_bcast( x_q,    ionode_id, intra_image_comm )

      ALLOCATE(lgamma_iq(nqs))
      DO iq=1,nqs
        lgamma_iq(iq)=(x_q(1,iq)==0.D0.AND.x_q(2,iq)==0.D0.AND.x_q(3,iq)==0.D0) 
      END DO

      IF (fpol) THEN
         CALL mp_bcast( nfs_,    ionode_id, intra_image_comm )
         IF (nfs_ /= nfs) &
            CALL errore('read_qu','wrong number of frequencies',1)

         CALL mp_bcast( fiu,    ionode_id, intra_image_comm )
      END IF

      RETURN
      !
    END SUBROUTINE read_qu

    SUBROUTINE read_partial_ph( irr, ierr )

    USE partial,    ONLY : done_irr
    USE efield_mod, ONLY : zstarue0_rec
    USE dynmat,     ONLY : dyn_rec
    USE control_ph, ONLY : trans, zue

    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: ierr
    INTEGER, INTENT(IN) :: irr

    !
    ierr=0
    IF (ionode) THEN
       IF (trans) THEN
          CALL iotk_scan_begin( iunpun, "PM_HEADER" )
          CALL iotk_scan_dat( iunpun,"DONE_IRR",done_irr(irr) )
          CALL iotk_scan_end( iunpun, "PM_HEADER" )
          CALL iotk_scan_begin( iunpun, "PARTIAL_MATRIX" )
          CALL iotk_scan_dat(iunpun,"PARTIAL_DYN", dyn_rec(:,:))
          IF (zue.AND.irr>0) CALL iotk_scan_dat(iunpun, &
                               "PARTIAL_ZUE", zstarue0_rec(:,:))

          CALL iotk_scan_end( iunpun, "PARTIAL_MATRIX" )
       ENDIF
    ENDIF
    IF (trans) THEN
       CALL mp_bcast( done_irr(irr),  ionode_id, intra_image_comm )
       CALL mp_bcast( dyn_rec,  ionode_id, intra_image_comm )
       IF (zue) CALL mp_bcast( zstarue0_rec,  ionode_id, intra_image_comm )
    ENDIF

    RETURN
    END SUBROUTINE read_partial_ph

    SUBROUTINE read_el_phon(irr, ierr)
    USE qpoint,     ONLY : nksq, nksqtot
    USE el_phon,    ONLY : el_ph_mat_rec, el_ph_mat_rec_col, done_elph, elph
    USE modes,      ONLY : npert
    USE wvfct,      ONLY : nbnd
    USE mp_pools,   ONLY : npool
    
    IMPLICIT NONE
    INTEGER,          INTENT(in) :: irr
    INTEGER,          INTENT(OUT) :: ierr
    REAL(DP) :: xkdum(3)
    INTEGER :: ik, npe, idum
    !
    ierr=0
    IF (.NOT. elph) RETURN
    npe=npert(irr)

    IF (npool>1) THEN
       ALLOCATE(el_ph_mat_rec_col(nbnd,nbnd,nksqtot,npe))
    ELSE
       el_ph_mat_rec_col => el_ph_mat_rec
    ENDIF

    IF (ionode) THEN
       CALL iotk_scan_begin(iunpun, "EL_PHON_HEADER")
       CALL iotk_scan_dat(iunpun, "DONE_ELPH", done_elph(irr))
       CALL iotk_scan_end(iunpun, "EL_PHON_HEADER")
       CALL iotk_scan_begin(iunpun, "PARTIAL_EL_PHON")
       CALL iotk_scan_dat(iunpun, "NUMBER_OF_K", idum)
       CALL iotk_scan_dat(iunpun, "NUMBER_OF_BANDS", idum)
       DO ik=1,nksqtot
          CALL iotk_scan_begin(iunpun, "K_POINT" // &
                                        TRIM( iotk_index( ik ) ) )
          CALL iotk_scan_dat(iunpun, "COORDINATES_XK", xkdum(:))
          CALL iotk_scan_dat(iunpun, "PARTIAL_ELPH", &
                                  el_ph_mat_rec_col(:,:,ik,:))
          CALL iotk_scan_end(iunpun, "K_POINT" // &
                                        TRIM( iotk_index( ik ) ) )
       ENDDO
       CALL iotk_scan_end(iunpun, "PARTIAL_EL_PHON")
    ENDIF

    CALL mp_bcast(done_elph(irr), ionode_id, intra_image_comm)
    CALL mp_bcast(el_ph_mat_rec_col, ionode_id, intra_image_comm)

    IF (npool > 1) THEN
       CALL el_ph_distribute(npe,el_ph_mat_rec,el_ph_mat_rec_col,&
                                                   nksqtot,nksq)
       DEALLOCATE(el_ph_mat_rec_col)
    ENDIF

    RETURN
    END SUBROUTINE read_el_phon

    SUBROUTINE read_modes( current_iq, ierr )
!
!   This routine reads the displacement patterns.
!
    USE modes, ONLY : nirr, npert, u, name_rap_mode, num_rap_mode
    USE el_phon, ONLY : elph 
    USE control_ph, ONLY : trans, zeu

    USE lr_symm_base, ONLY : nsymq, minus_q

    IMPLICIT NONE

    INTEGER,          INTENT(IN) :: current_iq
    INTEGER,          INTENT(OUT) :: ierr
    INTEGER :: imode0, imode, irr, ipert, iq, iu

    LOGICAL :: exst
    !
    ierr=0
    IF (ionode) THEN
       CALL iotk_scan_begin( iunpun, "IRREPS_INFO" )
       !
       CALL iotk_scan_dat(iunpun,"QPOINT_NUMBER",iq)
    ENDIF
    CALL mp_bcast( iq,  ionode_id, intra_image_comm )
    IF (iq /= current_iq) CALL errore('read_modes', &
              'problems with current_iq', 1 )

    IF (ionode) THEN

       CALL iotk_scan_dat(iunpun, "QPOINT_GROUP_RANK", nsymq)
       CALL iotk_scan_dat(iunpun, "MINUS_Q_SYM", minus_q)
       CALL iotk_scan_dat(iunpun, "NUMBER_IRR_REP", nirr)
       imode0=0
       DO irr=1,nirr
          CALL iotk_scan_begin( iunpun, "REPRESENTION"// &
                                     TRIM( iotk_index( irr ) ) )
          CALL iotk_scan_dat(iunpun,"NUMBER_OF_PERTURBATIONS", npert(irr))
          DO ipert=1,npert(irr)
             imode=imode0+ipert
             CALL iotk_scan_begin( iunpun, "PERTURBATION"// &
                                     TRIM( iotk_index( ipert ) ) )
             CALL iotk_scan_dat(iunpun,"SYMMETRY_TYPE_CODE", &
                                                        num_rap_mode(imode))
             CALL iotk_scan_dat(iunpun,"SYMMETRY_TYPE", name_rap_mode(imode))
             CALL iotk_scan_dat(iunpun,"DISPLACEMENT_PATTERN",u(:,imode))
             CALL iotk_scan_end( iunpun, "PERTURBATION"// &
                                     TRIM( iotk_index( ipert ) ) )
          ENDDO
          imode0=imode0+npert(irr)
          CALL iotk_scan_end( iunpun, "REPRESENTION"// &
                                     TRIM( iotk_index( irr ) ) )
       ENDDO
       !
       CALL iotk_scan_end( iunpun, "IRREPS_INFO" )
       !
    ENDIF

    CALL mp_bcast( nirr,  ionode_id, intra_image_comm )
    CALL mp_bcast( npert,  ionode_id, intra_image_comm )
    CALL mp_bcast( nsymq,  ionode_id, intra_image_comm )
    CALL mp_bcast( minus_q,  ionode_id, intra_image_comm )
    CALL mp_bcast( u,  ionode_id, intra_image_comm )
    CALL mp_bcast( name_rap_mode,  ionode_id, intra_image_comm )
    CALL mp_bcast( num_rap_mode,  ionode_id, intra_image_comm )

    RETURN
    END SUBROUTINE read_modes

    SUBROUTINE read_tensors( ierr )
!
!   This routine reads the tensors that have been already calculated 
!
    USE ions_base,  ONLY : nat
    USE control_ph, ONLY : done_epsil, done_start_zstar, done_zeu, done_zue
    USE ramanm,  ONLY : lraman, elop, ramtns, eloptns, done_lraman, done_elop
    USE efield_mod, ONLY : zstareu, zstarue, zstarue0, zstareu0, epsilon

    IMPLICIT NONE

    INTEGER,          INTENT(OUT) :: ierr
    INTEGER :: imode0, imode, ipol, irr, ipert, iq, iu, iunout
    !
    ierr=0
    IF (ionode) THEN
       CALL iotk_scan_begin( iunpun, "EF_TENSORS" )
       !
       CALL iotk_scan_dat( iunpun, "DONE_ELECTRIC_FIELD", done_epsil )
       CALL iotk_scan_dat( iunpun, "DONE_START_EFFECTIVE_CHARGE", done_start_zstar )
       CALL iotk_scan_dat( iunpun, "DONE_EFFECTIVE_CHARGE_EU", done_zeu )
       CALL iotk_scan_dat( iunpun, "DONE_EFFECTIVE_CHARGE_PH", done_zue )
       CALL iotk_scan_dat( iunpun, "DONE_RAMAN_TENSOR", done_lraman )
       CALL iotk_scan_dat( iunpun, "DONE_ELECTRO_OPTIC", done_elop )

       IF (done_epsil) &
          CALL iotk_scan_dat(iunpun,"DIELECTRIC_CONSTANT",epsilon)
       IF (done_start_zstar) &
          CALL iotk_scan_dat(iunpun,"START_EFFECTIVE_CHARGES",zstareu0)
       IF (done_zeu) &
          CALL iotk_scan_dat(iunpun,"EFFECTIVE_CHARGES_EU",zstareu)
       IF (done_lraman) &
          CALL iotk_scan_dat(iunpun,"RAMAN_TNS",ramtns)
       IF (done_elop) &
          CALL iotk_scan_dat(iunpun,"ELOP_TNS",eloptns)
       IF (done_zue) &
          CALL iotk_scan_dat(iunpun,"EFFECTIVE_CHARGES_UE",zstarue)
       !
       CALL iotk_scan_end( iunpun, "EF_TENSORS" )
       !
    ENDIF

    CALL mp_bcast( done_epsil,  ionode_id, intra_image_comm )
    CALL mp_bcast( done_start_zstar,  ionode_id, intra_image_comm )
    CALL mp_bcast( done_zeu,  ionode_id, intra_image_comm )
    CALL mp_bcast( done_zue,  ionode_id, intra_image_comm )
    CALL mp_bcast( done_lraman,  ionode_id, intra_image_comm )
    CALL mp_bcast( done_elop,  ionode_id, intra_image_comm )
    IF (done_epsil) CALL mp_bcast( epsilon, ionode_id, intra_image_comm )
    IF (done_start_zstar) THEN
       CALL mp_bcast( zstareu0, ionode_id, intra_image_comm )
       DO ipol=1,3
          DO imode=1,3*nat
             zstarue0(imode,ipol)=zstareu0(ipol,imode)
          ENDDO
       ENDDO
    ENDIF   
    IF (done_zeu) CALL mp_bcast( zstareu, ionode_id, intra_image_comm )
    IF (done_zue) CALL mp_bcast( zstarue, ionode_id, intra_image_comm )
    IF (done_lraman) CALL mp_bcast( ramtns, ionode_id, intra_image_comm )
    IF (done_elop) CALL mp_bcast( eloptns,  ionode_id, intra_image_comm )

    RETURN
    END SUBROUTINE read_tensors

  !----------------------------------------------------------------------------
    SUBROUTINE read_polarization( iu, ierr )
!
!   This routine reads the tensors that have been already calculated 
!
    USE ions_base,  ONLY : nat
    USE freq_ph, ONLY : fpol, done_iu, fiu, polar

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: iu
    INTEGER, INTENT(OUT) :: ierr
    !
    ierr=0
    IF ( .NOT.fpol ) RETURN
    IF (ionode) THEN
       CALL iotk_scan_begin( iunpun, "POLARIZ_IU" )
       !
       CALL iotk_scan_dat( iunpun, "FREQUENCY_IN_RY", fiu(iu) )
       CALL iotk_scan_dat( iunpun, "CALCULATED_FREQUENCY", &
                               done_iu(iu))
       IF (done_iu(iu)) CALL iotk_scan_dat( iunpun, &
                                     "POLARIZATION_IU", polar(:,:,iu) )
       !
       CALL iotk_scan_end( iunpun, "POLARIZ_IU" )
       !
    ENDIF

    CALL mp_bcast( fiu(iu),  ionode_id, intra_image_comm )
    CALL mp_bcast( done_iu(iu),  ionode_id, intra_image_comm )
    IF ( done_iu(iu) ) &
             CALL mp_bcast( polar(:,:,iu),  ionode_id, intra_image_comm )

    RETURN
    END SUBROUTINE read_polarization

  !----------------------------------------------------------------------------
    SUBROUTINE check_directory_phsave(  )
  !----------------------------------------------------------------------------
  ! ...
  ! ... This routine sets the situation of the grid according to
  ! ... the files that it finds on the directory .phsave.
  ! ... Check if representation files exist and which representations 
  ! ... have been already calculated.
  ! ... set the initial information on the grid
  ! ... it sets done_irr_iq to .true. for the q and the 
  ! ... representations that have already been done.
  ! ... Moreover it sets irr_iq, the number of representations for each q,
  ! ... nsymq_iq the size of the small group of each q and npert_irr_iq
  ! ... the number of perturbations for each irr and q.
  !
  USE kinds, ONLY : DP
  USE disp, ONLY : nqs, done_iq
  USE grid_irr_iq, ONLY : comp_irr_iq, done_irr_iq, irr_iq, done_elph_iq
  USE control_ph, ONLY : trans, current_iq, low_directory_check
  USE el_phon,    ONLY : elph
  ! 
  IMPLICIT NONE
  !
  CHARACTER(LEN=256)  :: dirname, filename, filename1
  INTEGER :: iunout, iq, irr, ierr
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  !
  ierr=0
  IF (ionode) THEN
     CALL iotk_free_unit( iunout, ierr )
  ENDIF
  !
  CALL mp_bcast( ierr, ionode_id, intra_image_comm )
  !
  CALL errore( 'check_directory_phsave', &
                   'no free units to write or read ', ierr )

  dirname = TRIM( tmp_dir_ph ) // TRIM( prefix ) // '.phsave'
  ierr=0
  DO iq=1, nqs
     IF ( ionode ) THEN
        IF (trans.OR.elph) THEN
!
!    NB: the representation 0 is the initial dynamical matrix calculated by 
!        dyn0. If it finds the file read the relevant information
!
           filename= TRIM( dirname ) // '/dynmat.' // &
                     TRIM(int_to_char(iq)) // '.' 

           DO irr=0,irr_iq(iq)
              IF (comp_irr_iq(irr,iq).OR..NOT.low_directory_check) THEN
                 filename1=TRIM(filename) // TRIM(int_to_char(irr)) // '.xml'
                 INQUIRE(FILE=TRIM(filename1), EXIST=exst)
                 IF (.NOT.exst) CYCLE
                 CALL iotk_open_read(iunout, FILE = TRIM(filename1), &
                                       BINARY = .FALSE., IERR = ierr )
                 IF (ierr /= 0 ) GOTO 100
                 CALL iotk_scan_begin( iunout, "PM_HEADER" )
                 CALL iotk_scan_dat(iunout,"DONE_IRR",done_irr_iq(irr,iq))
                 CALL iotk_scan_end( iunout, "PM_HEADER" )
                 CALL iotk_close_read(iunout)
              ENDIF
           END DO
!
!   Check for the electron phonon files
!
           IF (elph) THEN
              filename= TRIM( dirname ) // '/elph.' // &
                        TRIM(int_to_char(iq)) // '.' 

              DO irr=1,irr_iq(iq)
                 IF (comp_irr_iq(irr,iq).OR..NOT.low_directory_check) THEN
                    filename1=TRIM(filename) // TRIM(int_to_char(irr)) // '.xml'
                    INQUIRE(FILE=TRIM(filename1), EXIST=exst)
                    IF (.NOT.exst) CYCLE
                    CALL iotk_open_read(iunout, FILE = TRIM(filename1), &
                                          BINARY = .FALSE., IERR = ierr )
                    IF (ierr /= 0 ) GOTO 100
                    CALL iotk_scan_begin(iunout, "EL_PHON_HEADER")
                    CALL iotk_scan_dat(iunout, "DONE_ELPH", done_elph_iq(irr,iq))
                    CALL iotk_scan_end(iunout, "EL_PHON_HEADER")
                    CALL iotk_close_read(iunout)
                 ENDIF
              ENDDO
           END IF
        END IF 
        done_iq(iq)=.TRUE.
        DO irr=1,irr_iq(iq)
           IF (comp_irr_iq(irr,iq).AND..NOT.done_irr_iq(irr,iq)) &
                                        done_iq(iq)=.FALSE.
           IF (elph) THEN
              IF (comp_irr_iq(irr,iq).AND..NOT.done_elph_iq(irr,iq)) &
                                        done_iq(iq)=.FALSE.
           ENDIF
        ENDDO
        IF (comp_irr_iq(0,iq).AND..NOT.done_irr_iq(0,iq)) done_iq(iq)=.FALSE.
     END IF
  END DO
100 CALL mp_bcast( ierr, ionode_id, intra_image_comm )
  IF (ierr /= 0) CALL errore('check_directory_phsave','opening file',1)
  !
  CALL mp_bcast( done_iq, ionode_id, intra_image_comm )
  CALL mp_bcast( done_irr_iq, ionode_id, intra_image_comm )
  IF (elph) CALL mp_bcast( done_elph_iq, ionode_id, intra_image_comm )
  !
  RETURN
  !
  END SUBROUTINE check_directory_phsave

  !----------------------------------------------------------------------------
    SUBROUTINE check_available_bands(  )
  !----------------------------------------------------------------------------
  ! ...
  ! ... This routine checks which bands are available on disk and
  ! ... sets the array done_bands(iq) to .true. for each q point
  ! ... for which the bands are present.
  ! ... If lqdir is .false. only the bands corresponding to current_iq
  ! ... can be present, whereas if lqdir is .true. several q points
  ! ... might have calculated the bands and saved them on disk.
  !
  USE kinds, ONLY : DP
  USE disp, ONLY : nqs, x_q, lgamma_iq
  USE io_files, ONLY : tmp_dir
  USE control_ph, ONLY : tmp_dir_ph, lqdir, current_iq, newgrid
  USE grid_irr_iq, ONLY : done_bands
  ! 
  IMPLICIT NONE
  !
  CHARACTER(LEN=256)  :: dirname, filename, dir_phq, tmp_dir_save
  INTEGER :: iq
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  LOGICAL :: lgamma, exst, exst_restart, exst_recover
  !
  ! We check if the file data-file.xml is present 
  ! in the directory where it should be. If lqdir=.false. only the bands 
  ! of current_iq might be present, otherwise we have to check all q points.
  ! If the file is present and there is a restart file, the bands are not
  ! done yet.
  ! For the gamma point done_bands might be false only with newgrid. 
  !
  done_bands=.FALSE.
  dirname = TRIM( tmp_dir_ph ) // TRIM( prefix ) // '.save'
  tmp_dir_save=tmp_dir

  DO iq=1, nqs
     IF ( iq == current_iq .OR. lqdir) THEN
        IF (lqdir .AND. .NOT. lgamma_iq(iq)) THEN
           dir_phq= TRIM (tmp_dir_ph) // TRIM(prefix) //  &
                       & '.q_' // TRIM(int_to_char(iq)) // '/'
           dirname= TRIM (dir_phq) &
                       &  //TRIM(prefix)//'.save'
           tmp_dir=dir_phq
        ELSE
           tmp_dir=tmp_dir_ph
        ENDIF
        !
        filename=TRIM(dirname) // '/data-file.xml'
        !
        IF (ionode) inquire (file =TRIM(filename), exist = exst)
        !
        CALL mp_bcast( exst, ionode_id, intra_image_comm )
        !
        exst_restart=.FALSE.
        IF (exst) CALL check_restart_recover(exst_recover, exst_restart)
        !
        IF (exst.AND..NOT.exst_restart) done_bands(iq)=.TRUE.
     END IF
     IF (lgamma_iq(iq).AND..NOT.newgrid) done_bands(iq) = .TRUE.
  END DO
  tmp_dir=tmp_dir_save
  !
  RETURN
  !
  END SUBROUTINE check_available_bands

   SUBROUTINE allocate_grid_variables()
!
!  This routine allocates and initializes the grid variables when the 
!  nqs and x_q have been decided, either reading them from file when
!  recover is .true. or recalculating them from scratch  
!
   USE disp, ONLY : nqs, done_iq, comp_iq, omega_disp              
   USE grid_irr_iq, ONLY : done_irr_iq, irr_iq, nsymq_iq, &
                           comp_irr_iq, npert_irr_iq, done_bands, &
                           done_elph_iq
   USE freq_ph, ONLY : done_iu, comp_iu, nfs
   USE ions_base, ONLY : nat
   USE el_phon, ONLY : elph_simple, gamma_disp, el_ph_nsigma
   USE control_ph, ONLY : qplot

   IMPLICIT NONE

   ALLOCATE(done_iq(nqs))
   ALLOCATE(done_bands(nqs))
   ALLOCATE(comp_iq(nqs))
   ALLOCATE(irr_iq(nqs))
   ALLOCATE(done_irr_iq(0:3*nat,nqs))
   ALLOCATE(done_elph_iq(1:3*nat,nqs))
   ALLOCATE(comp_irr_iq(0:3*nat,nqs))
   ALLOCATE(nsymq_iq(nqs))
   ALLOCATE(npert_irr_iq(3*nat,nqs))

   ALLOCATE(done_iu(nfs))
   ALLOCATE(comp_iu(nfs))

   done_iq=.FALSE.
   done_bands=.FALSE.
   done_irr_iq=.FALSE.
   done_elph_iq=.FALSE.
   done_iu=.FALSE.

   comp_iu=.TRUE.
   comp_iq=.TRUE.
   comp_irr_iq=.TRUE.

   irr_iq=3*nat
   nsymq_iq=0
   npert_irr_iq=0

   IF (qplot) THEN
      ALLOCATE(omega_disp(3*nat,nqs))
      IF (elph_simple) ALLOCATE(gamma_disp(3*nat,el_ph_nsigma,nqs))
   ENDIF

   RETURN
   END SUBROUTINE allocate_grid_variables

   SUBROUTINE destroy_status_run()
   USE start_k,     ONLY : xk_start, wk_start
   USE disp,        ONLY : nqs, x_q, done_iq, comp_iq, lgamma_iq, omega_disp
   USE grid_irr_iq, ONLY : done_irr_iq, irr_iq, nsymq_iq, &
                          npert_irr_iq, comp_irr_iq, done_bands, done_elph_iq
   USE el_phon,     ONLY : gamma_disp
   USE freq_ph,     ONLY : comp_iu, done_iu, fiu
   IMPLICIT NONE

   IF (ALLOCATED(x_q)) DEALLOCATE(x_q)
   IF (ALLOCATED(lgamma_iq)) DEALLOCATE(lgamma_iq)
   IF (ALLOCATED(done_bands)) DEALLOCATE(done_bands)
   IF (ALLOCATED(done_iq)) DEALLOCATE(done_iq)
   IF (ALLOCATED(comp_iq)) DEALLOCATE(comp_iq)
   IF (ALLOCATED(irr_iq)) DEALLOCATE(irr_iq)
   IF (ALLOCATED(done_irr_iq)) DEALLOCATE(done_irr_iq)
   IF (ALLOCATED(done_elph_iq)) DEALLOCATE(done_elph_iq)
   IF (ALLOCATED(comp_irr_iq)) DEALLOCATE(comp_irr_iq)
   IF (ALLOCATED(nsymq_iq)) DEALLOCATE(nsymq_iq)
   IF (ALLOCATED(npert_irr_iq)) DEALLOCATE(npert_irr_iq)
   IF (ALLOCATED(fiu)) DEALLOCATE(fiu)
   IF (ALLOCATED(done_iu)) DEALLOCATE(done_iu)
   IF (ALLOCATED(comp_iu)) DEALLOCATE(comp_iu)
   IF (ALLOCATED(omega_disp)) DEALLOCATE(omega_disp)
   IF (ALLOCATED(gamma_disp)) DEALLOCATE(gamma_disp)
!
! Note that these two variables are allocated by read_file. 
! They cannot be deallocated by clean_pw because the starting xk and wk 
! points must be known at each q point.
! The logic of these two variables must be improved.
!
  IF (ALLOCATED( xk_start )) DEALLOCATE( xk_start )
  IF (ALLOCATED( wk_start )) DEALLOCATE( wk_start )

   END SUBROUTINE destroy_status_run

   SUBROUTINE ph_restart_set_filename( what, irr, current_iq, iflag, ierr)
!
!    This subroutine sets the filename for each action required by what
!    and opens the appropriate file for reading or writing
!
      USE io_global,    ONLY : ionode, ionode_id
      USE xml_io_base,  ONLY : create_directory
      USE freq_ph,      ONLY : fpol
      USE mp_images,    ONLY : intra_image_comm
      USE mp,           ONLY : mp_bcast

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: irr, current_iq, iflag
      INTEGER, INTENT(OUT) :: ierr
      CHARACTER(LEN=*), INTENT(IN) :: what

      CHARACTER(LEN=256)  :: dirname, filename
      CHARACTER(LEN=6) :: int_to_char
      LOGICAL :: exst

      ierr=0
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
      CALL errore( 'ph_restart_set_filename', &
                   'no free units to write or read ', ierr )
      !
      dirname = TRIM( tmp_dir_ph ) // TRIM( prefix ) // '.phsave'
      !
      ! ... create the main restart directory
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
            filename = TRIM( dirname ) // '/control_ph.xml'
            IF (iflag==1) THEN
               CALL iotk_open_write( iunpun, FILE = TRIM(filename), &
                                             BINARY = .FALSE., IERR = ierr )
            ELSE
               INQUIRE( FILE=TRIM(filename), EXIST=exst )
               IF (.NOT.exst) GOTO 100
               CALL iotk_open_read( iunpun, FILE = TRIM(filename), &
                                            BINARY = .FALSE., IERR = ierr )
            ENDIF
         ELSEIF (what=='status_ph') THEN
            filename=TRIM( dirname ) //'/status_run.xml'
            IF (iflag==1) THEN
               CALL iotk_open_write( iunpun, FILE = TRIM( filename ),  &
                                            BINARY = .FALSE., IERR = ierr )
            ELSE
               INQUIRE( FILE=TRIM(filename), EXIST=exst )
               IF (.NOT.exst) GOTO 100
               CALL iotk_open_read( iunpun, FILE = TRIM( filename ),  &
                                            BINARY = .FALSE., IERR = ierr )
            ENDIF
         ELSEIF (what=='data_u') THEN
            filename= TRIM( dirname ) // '/patterns.' // &
                      TRIM(int_to_char(current_iq)) // '.xml'
            IF (iflag==1) THEN
               CALL iotk_open_write( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
            ELSE
               INQUIRE( FILE=TRIM(filename), EXIST=exst )
               IF (.NOT.exst) GOTO 100
               CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
            ENDIF
         ELSEIF (what=='data_dyn') THEN
            filename= TRIM( dirname ) // '/dynmat.' // &
                      TRIM(int_to_char(current_iq)) // '.' //  &
                      TRIM(int_to_char(irr)) // '.xml'
            IF (iflag==1) THEN
               CALL iotk_open_write( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
            ELSE
               INQUIRE( FILE=TRIM(filename), EXIST=exst )
               IF (.NOT.exst) GOTO 100
               CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
            ENDIF
         ELSEIF (what=='tensors') THEN
            filename= TRIM( dirname ) // '/tensors.xml'
            IF (iflag==1) THEN
               CALL iotk_open_write( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
            ELSE
               INQUIRE( FILE=TRIM(filename), EXIST=exst )
               IF (.NOT.exst) GOTO 100
               CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
            ENDIF
         ELSEIF (what=='polarization') THEN
            IF (.NOT. fpol) RETURN
            filename= TRIM( dirname ) // '/polarization.'// &
                      TRIM(int_to_char(irr)) // '.xml'
            IF (iflag==1) THEN
               CALL iotk_open_write( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
            ELSE
               INQUIRE( FILE=TRIM(filename), EXIST=exst )
               IF (.NOT.exst) GOTO 100
               CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
            ENDIF
         ELSEIF (what=='el_phon') THEN
            filename= TRIM( dirname ) // '/elph.' // &
                      TRIM(int_to_char(current_iq)) // '.' //  &
                      TRIM(int_to_char(irr)) // '.xml'
            IF (iflag==1) THEN
               CALL iotk_open_write( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
            ELSE
               INQUIRE( FILE=TRIM(filename), EXIST=exst )
               IF (.NOT.exst) GOTO 100
               CALL iotk_open_read( iunpun, FILE = TRIM( filename ), &
                                          BINARY = .FALSE., IERR = ierr )
            ENDIF
         ELSE
            ierr = 1
         ENDIF
         !
      END IF
100   IF (iflag /= 0) THEN
         CALL mp_bcast( exst, ionode_id, intra_image_comm )
!
!     If the file does not exist and we must read from it, we return with
!     an error message.
!
         IF (.NOT.exst) THEN
            ierr=100
            RETURN
         ENDIF
      ENDIF
      !
      CALL mp_bcast( ierr, ionode_id, intra_image_comm )
      !
      CALL errore( 'ph_restart_set_filename ', &
                   'cannot open file for reading or writing', ierr )
  
    RETURN
    END SUBROUTINE ph_restart_set_filename
    !
END MODULE ph_restart


