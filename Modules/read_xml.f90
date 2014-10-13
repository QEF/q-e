!
!---------------------------------------------------------!
!   This module handles the reading of fields and cards   !
!   in case of xml input                                  !
!                                                         !
!   written by Simone Ziraldo (08/2010)                   !
!---------------------------------------------------------!
MODULE read_xml_module
  !
  !
  USE input_parameters
  !
  USE io_global, ONLY : ionode, ionode_id, xmlinputunit => qestdin
  USE mp,        ONLY : mp_bcast
  USE mp_images, ONLY : intra_image_comm
  USE iotk_module, ONLY : iotk_attlenx
  !
  ! ...default and checkin of fields
  !
  USE read_namelists_module, ONLY : control_defaults, system_defaults,&
       electrons_defaults, wannier_ac_defaults, ions_defaults, &
       cell_defaults, press_ai_defaults, wannier_defaults, control_bcast, &
       system_bcast, electrons_bcast, ions_bcast,cell_bcast, &
       press_ai_bcast, wannier_bcast, wannier_ac_bcast, control_checkin, &
       system_checkin, electrons_checkin, ions_checkin, cell_checkin, &
       wannier_checkin, wannier_ac_checkin, fixval
  !
  !
  USE read_xml_fields_module, ONLY : read_xml_fields
  USE read_xml_cards_module, ONLY : card_xml_atomic_species, card_xml_atomic_list, &
       card_xml_cell, card_xml_kpoints, card_xml_occupations, &
       card_xml_constraints, card_xml_plot_wannier, card_default, card_bcast
  !
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: read_xml
  !
CONTAINS
  !
  !
  !--------------------------------------------------------!
  !    This routine organizes the reading of the xml file  !
  !    depending on the program                            !
  !--------------------------------------------------------!
  SUBROUTINE read_xml( prog, attr )
    !
    !
    IMPLICIT NONE
    !
    !
    CHARACTER(len = 2), INTENT(IN) :: prog
    CHARACTER(len = *), INTENT(IN) :: attr
    INTEGER :: ierr
    !
    SELECT CASE (prog)
       !
    CASE ('PW')
       !
       CALL read_xml_common( attr, 'PW' )
       CALL read_xml_pw()
       !
    CASE ('CP')
       !
       CALL read_xml_common( attr, 'CP' )
       CALL read_xml_cp()
       !
    CASE default
       !
       CALL errore('read_xml', "xml input isn't implemented for "//prog//' program', 1)
       !
    END SELECT
    !
    !
    RETURN
    !
  END SUBROUTINE read_xml
  !
  !
  !--------------------------------------------------------!
  ! Common part of the reading: setting default values,    !
  ! reading of cell and atomic_species cards               !
  !--------------------------------------------------------!
  SUBROUTINE read_xml_common( attr, prog )
    !
    !
    USE iotk_module, ONLY : iotk_scan_attr
    !
    !
    IMPLICIT NONE
    !
    !
    CHARACTER (len = *), INTENT(IN) :: attr, prog
    !
    CHARACTER (len = 256) :: dummy
    INTEGER :: ierr
    LOGICAL :: found
    !
    !
    ! ... default settings for all parameters
    !
    CALL control_defaults( prog )
    CALL system_defaults( prog )
    CALL electrons_defaults( prog )
    CALL ions_defaults( prog )
    CALL cell_defaults( prog )
    CALL wannier_defaults( prog )
    CALL wannier_ac_defaults( prog )
    !
    !
    ! ... reading the attributes of the xml root node
    !
    IF (ionode) THEN
       !
       CALL iotk_scan_attr( attr, 'calculation', dummy, found = found, ierr = ierr )
       IF ( .not. found ) CALL errore( 'read_xml_common', 'attribute calculation of root &
            &node is compulsory', abs(ierr) )
       !
       IF ( ierr /= 0 ) CALL errore( 'read_xml_common', 'error reading calculation &
            &attribute of root node', 1 )
       calculation = trim( dummy )
       !
       CALL iotk_scan_attr( attr, 'prefix', dummy, found = found, ierr = ierr )
       IF ( ierr /= 0 ) CALL errore( 'read_xml_common', 'error reading prefix attribute &
            &of root node', abs(ierr) )
       IF ( found ) prefix = trim( dummy )
       !
       CALL iotk_scan_attr( attr, 'title', dummy, found = found, ierr = ierr )
       IF ( ierr /= 0 ) CALL errore( 'read_xml_common', 'error reading title attribute &
            &of root node', 1 )
       IF ( found ) title = trim( dummy )
       !
    END IF
    !
    !  ... bcast the read attributes
    !
    CALL mp_bcast( calculation, ionode_id, intra_image_comm )
    CALL mp_bcast( prefix, ionode_id, intra_image_comm )
    CALL mp_bcast( title, ionode_id, intra_image_comm )
    
    !
    ! ... fixing some default values using the calculation variable
    !
    CALL fixval( prog )
    !
    ! ... why this is compulsory? ( read autopilot.f90 )
    CALL card_default( 'INIT_AUTOPILOT' )
    !
    !
    !   ... reading CELL card
    !   
    CALL card_default( 'CELL' )
    !
    IF ( ionode ) THEN
       !
       CALL card_xml_cell( )
       !
    END IF
    !
    CALL card_bcast( 'CELL' )
    !
    !
    ! ...  reading ATOMIC_SPECIES card
    !
    CALL card_default( 'ATOMIC_SPECIES' )
    !
    IF ( ionode ) THEN
       !
       CALL card_xml_atomic_species( )
       !
    END IF
    !
    CALL card_bcast( 'ATOMIC_SPECIES' )
    !
    RETURN
    !
  END SUBROUTINE read_xml_common
  !
  !
  !--------------------------------------------------------!
  ! The remaining part of the reading for PW: fields and   !
  ! other cards                                            !
  !--------------------------------------------------------!
  SUBROUTINE read_xml_pw( )
    !
    !
    USE iotk_module, ONLY : iotk_scan_begin, iotk_scan_end
    USE iotk_unit_interf, ONLY : iotk_rewind
    !
    !
    IMPLICIT NONE
    !
    !
    INTEGER :: ierr
    CHARACTER (len = iotk_attlenx) :: attr
    CHARACTER (len = 30) :: field, card
    LOGICAL :: found_al, found
    !
    !
    ! ... reading ATOMIC_LIST card
    !
    CALL card_default( 'ATOMIC_LIST' )
    !
    IF ( ionode ) THEN
       !
       CALL iotk_scan_begin( xmlinputunit, 'atomic_list', found = found_al, ierr = ierr )
       IF ( ierr /= 0 ) CALL errore( 'read_xml_pw', 'error scanning begin &
            &of atomic_list card', abs(ierr) )
       !
       IF ( found_al ) THEN
          !
          CALL iotk_scan_end( xmlinputunit, 'atomic_list', ierr = ierr )
          IF ( ierr /= 0 ) CALL errore( 'read_xml_pw', 'error scanning end &
               &of atomic_list card', abs( ierr ) )
          !
          CALL card_xml_atomic_list( )
          !
       ELSE
          !
          CALL errore('read_xml_pw',"card atomic_list is missing", 1 )
          !
       ENDIF
    ENDIF
    !
    CALL mp_bcast( found_al, ionode_id, intra_image_comm)
    !
    CALL card_bcast( 'ATOMIC_LIST' )
    !
    ! ... reading all the FIELDS
    !
    !
    ! ... we need to know if startingwfc and starting pot are set
    startingwfc = 'none'
    startingpot = 'none'
    !
    IF (ionode) THEN
       !
       CALL read_xml_fields()
       !
    END IF
    !
    !
    ! ... some fixval that the previous call of fixval wasn't
    ! ... able to do
    !
    IF ( calculation == 'nscf' .or. calculation == 'bands'  ) THEN
       !
       IF (startingpot == 'none') startingpot = 'file'
       IF (startingwfc == 'none') startingwfc = 'atomic+random'
       !
    ELSE IF ( restart_mode == 'from_scratch' ) THEN
       !
       IF (startingwfc == 'none') startingwfc = 'atomic+random'
       IF (startingpot == 'none') startingpot = 'atomic'
       !
    ELSE
       !
       IF (startingwfc == 'none') startingwfc = 'file'
       IF (startingpot == 'none') startingpot = 'file'
       !
    END IF
    !
    !
    ! 
    ! ... checkin of all the parameters inserted in the fields
    !
    IF ( ionode ) THEN
       !
       CALL control_checkin( 'PW' )
       CALL system_checkin( 'PW' )
       CALL electrons_checkin( 'PW' )
       CALL ions_checkin( 'PW' )
       CALL cell_checkin( 'PW' )
       CALL wannier_checkin( 'PW' )
       CALL wannier_ac_checkin( 'PW' )
       !
    END IF
    !
    !
    ! ... bcast all the field parameters
    !
    CALL control_bcast( )
    CALL system_bcast( )
    CALL electrons_bcast( )
    CALL ions_bcast( )
    CALL cell_bcast()
    CALL press_ai_bcast()
    CALL wannier_bcast()
    CALL wannier_ac_bcast()
    !
    !
    ! ... second step : reading of the remaining cards
    !
    !
    ! ... reading CONSTRAINTS card
    !
    card = 'constraints'
    CALL card_default( 'CONSTRAINTS' )
    !
    IF ( ionode ) THEN
       !
       CALL iotk_scan_begin( xmlinputunit, trim(card), found = found, ierr = ierr )
       IF ( ierr /= 0 ) GO TO 9
       !
       IF ( found ) THEN
          !
          CALL card_xml_constraints( )
          !
          CALL iotk_scan_end( xmlinputunit, trim(card), ierr = ierr)
          IF ( ierr /= 0 ) GOTO 10
          !
       ELSE
          !
          ! ... due to a iotk problem with gfortran compiler
          CALL iotk_rewind( xmlinputunit )
          !
       END IF
       !
    END IF
    !
    CALL mp_bcast ( found, ionode_id, intra_image_comm )
    !
    IF ( found ) CALL card_bcast( 'CONSTRAINTS' )
    !
    !
    ! ... reading K_POINTS card
    !
    card = 'k_points'
    CALL card_default( 'K_POINTS' )
    !
    IF ( ionode ) THEN
       !
       CALL iotk_scan_begin( xmlinputunit, trim( card ), attr = attr, found = found,&
            ierr = ierr )
       IF ( ierr /= 0 ) GO TO 9
       !
       IF ( found ) THEN
          !
          CALL card_xml_kpoints( attr )
          !
          CALL iotk_scan_end( xmlinputunit, trim( card ), ierr = ierr)
          IF ( ierr /= 0 ) GOTO 10
          !
       ELSE
          !
          CALL errore('read_xml_pw', 'K_POINTS card was not found', 1)
          !
       END IF
       !
    END IF
    !
    CALL card_bcast( 'K_POINTS' )
    !
    !
    ! ... reading OCCUPATIONS card
    !
    card = 'occupations'
    CALL card_default( 'OCCUPATIONS' )
    !
    IF ( ionode ) THEN
       !
       CALL iotk_scan_begin( xmlinputunit, trim( card ), found = found, ierr = ierr )
       IF ( ierr /= 0 ) GO TO 9
       !
       IF  ( found ) THEN
          !
          CALL card_xml_occupations()
          !
          CALL iotk_scan_end( xmlinputunit, trim( card ), ierr = ierr )
          IF ( ierr /= 0 ) GOTO 10
          !
       ELSE
          !
          ! ... due to a iotk problem with gfortran compiler
          CALL iotk_rewind( xmlinputunit )
          !
       END IF
       !
    END IF
    !
    CALL mp_bcast ( found, ionode_id, intra_image_comm )
    !
    IF ( found ) CALL card_bcast( 'OCCUPATIONS' )
    !
    RETURN
    !
9   CALL errore('read_xml_pw', 'error reading begin tag of '//card//' card', ABS( ierr ) )
10  CALL errore('read_xml_pw', 'error reading end tag of '//card//' card', ABS( ierr ) )
    !
    !
  END SUBROUTINE read_xml_pw
  !
  !
  !
  !--------------------------------------------------------!
  ! The rest of the reading for CP program : fileds and    !
  ! other cards                                            !
  !--------------------------------------------------------!
  SUBROUTINE read_xml_cp( )
    !
    !
    USE iotk_module, ONLY : iotk_scan_begin, iotk_scan_end
    USE iotk_unit_interf, ONLY : iotk_rewind
    !
    !
    IMPLICIT NONE
    !
    !
    INTEGER :: ierr
    CHARACTER (len = iotk_attlenx) :: attr
    CHARACTER (len = 30) :: field, card
    LOGICAL :: found
    !
    !
    ! ... reading ATOMIC_LIST cards
    !
    !
    CALL card_default( 'ATOMIC_LIST' )
    !
    IF ( ionode ) THEN
       !
       CALL card_xml_atomic_list ( )
       !
    END IF
    !
    CALL card_bcast( 'ATOMIC_LIST' )
    !
    !
    ! ... reading all the FIELDS
    !
    IF (ionode) THEN
       !
       CALL read_xml_fields()
       !
    END IF
    !
    ! 
    ! ... checkin of all the parameters inserted in the fields
    !
    IF ( ionode ) THEN
       !
       CALL control_checkin( 'CP' )
       CALL system_checkin( 'CP' )
       CALL electrons_checkin( 'CP' )
       CALL ions_checkin( 'CP' )
       CALL cell_checkin( 'CP' )
       CALL wannier_checkin( 'CP' )
       CALL wannier_ac_checkin( 'CP' )
       !
    END IF
    !
    !
    ! ... bcast all the field parameters
    !
    CALL control_bcast( )
    CALL system_bcast( )
    CALL electrons_bcast( )
    CALL ions_bcast( )
    CALL cell_bcast()
    CALL press_ai_bcast()
    CALL wannier_bcast()
    CALL wannier_ac_bcast()
    !
    !
    ! ... second step : reading of the remaining cards
    !
    !
    ! ... reading CONSTRAINTS card
    !
    card = 'constraints'
    CALL card_default( 'CONSTRAINTS' )
    !
    IF ( ionode ) THEN
       !
       CALL iotk_scan_begin( xmlinputunit, trim(card), found = found, ierr = ierr )
       IF ( ierr /= 0 ) GO TO 9
       !
       IF ( found ) THEN
          !
          CALL card_xml_constraints( )
          !
          CALL iotk_scan_end( xmlinputunit, trim(card), ierr = ierr)
          IF ( ierr /= 0 ) GOTO 10
          !
       ELSE
          !
          ! ... due to a iotk problem with gfortran compiler
          CALL iotk_rewind( xmlinputunit )
          !
       END IF
       !
    END IF
    !
    CALL mp_bcast ( found, ionode_id, intra_image_comm )
    !
    IF ( found ) CALL card_bcast( 'CONSTRAINTS' )
    !
    ! ... reading OCCUPATIONS card
    !
    card = 'occupations'
    CALL card_default( 'OCCUPATIONS' )
    !
    IF ( ionode ) THEN
       !
       CALL iotk_scan_begin( xmlinputunit, trim( card ), found = found, ierr = ierr )
       IF ( ierr /= 0 ) GO TO 9
       !
       IF  ( found ) THEN
          !
          CALL card_xml_occupations()
          !
          CALL iotk_scan_end( xmlinputunit, trim( card ), ierr = ierr )
          IF ( ierr /= 0 ) GOTO 10
          !
       ELSE
          !
          ! ... due to a iotk problem with gfortran compiler
          CALL iotk_rewind( xmlinputunit )
          !
       END IF
       !
    END IF
    !
    CALL mp_bcast ( found, ionode_id, intra_image_comm )
    !
    IF ( found ) CALL card_bcast( 'OCCUPATIONS' )
    !
    card = 'plot_wannier'
    CALL card_default( 'PLOT_WANNIER' )
    !
    IF ( ionode ) THEN
       !
       CALL iotk_scan_begin( xmlinputunit, trim( card ), found = found, ierr = ierr )
       IF ( ierr /= 0 ) GO TO 9
       !
       IF ( found ) THEN
          !
          CALL card_xml_plot_wannier()
          !
          CALL iotk_scan_end( xmlinputunit, trim( card ), ierr = ierr )
          IF ( ierr /= 0 ) GOTO 10
          !
       ELSE
          !
          ! ... due to a iotk problem with gfortran compiler
          CALL iotk_rewind( xmlinputunit )
          !
       END IF
       !
    END IF
    !
    CALL mp_bcast ( found, ionode_id, intra_image_comm )
    !
    IF ( found ) CALL card_bcast( 'PLOT_WANNIER' )
    !
    !
    !
    !
    RETURN
    !
9   CALL errore('read_xml_cp', 'error reading begin tag of '//card//' card', ABS( ierr ) )
10  CALL errore('read_xml_cp', 'error reading end tag of '//card//' card', ABS( ierr ) )
    !
    !
  END SUBROUTINE read_xml_cp
  !
  !
  !
END MODULE read_xml_module
