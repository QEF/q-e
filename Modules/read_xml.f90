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
  USE io_files, ONLY : xmlinputunit
  !
  ! ...default and checkin of fields
  !
  USE read_namelists_module, ONLY : control_defaults, system_defaults,&
       ee_defaults, electrons_defaults, wannier_ac_defaults, ions_defaults, &
       cell_defaults, press_ai_defaults, wannier_defaults, control_bcast, &
       system_bcast, ee_bcast, electrons_bcast, ions_bcast,cell_bcast, &
       press_ai_bcast, wannier_bcast, wannier_ac_bcast, control_checkin, &
       system_checkin, electrons_checkin, ions_checkin, cell_checkin, &
       wannier_checkin, wannier_ac_checkin, fixval
  !
  !
  USE read_xml_fields_module, ONLY : read_xml_fields
  USE read_xml_cards_module, ONLY : card_xml_atomic_species, card_xml_atomic_list, &
       card_xml_chain, card_xml_cell, card_xml_kpoints, card_xml_occupations, &
       card_xml_constraints, card_xml_climbing_images, card_default, card_bcast
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
  !    This routine sets default values of the parameters  !
  !    and organizes the reading of xml input              !
  !--------------------------------------------------------!
  SUBROUTINE read_xml( prog, step, attr )
    !
    !
    IMPLICIT NONE
    !
    !
    CHARACTER(len = 2), INTENT(IN) :: prog
    CHARACTER(len = *), INTENT(IN), OPTIONAL :: attr
    INTEGER, INTENT(IN) :: step
    !
    !
    ! ... default settings for all parameters
    !
    IF (step == 1) THEN
       CALL control_defaults( prog )
       CALL system_defaults( prog )
       CALL electrons_defaults( prog )
       CALL ions_defaults( prog )
       CALL cell_defaults( prog )
       CALL ee_defaults( prog )
       CALL wannier_defaults( prog )
       CALL wannier_ac_defaults( prog )
    END IF
    !
    ! ... begins the reading of xml input file
    !
    SELECT CASE (prog)
       !
    CASE ('PW')
       !
       CALL read_xml_pw( attr, step )
       !
    CASE default
       CALL errore('read_xml', "xml input isn't implemented for "//prog//' program', 1)
       !
       !
    END SELECT
    !
    RETURN
    !
  END SUBROUTINE read_xml
  !
  !
  !
  !
  !------------------------------------------------------------------------!
  !    This routine reads the xml file for the PW code                     !
  !    it has two steps:                                                   !
  !    1) reading of CELL, ATOMIC_SPECIES and ATOMIC_POSITIONS cells       !
  !       and all the fields                                               !
  !    2) reading of the remaining cards                                   !
  !    this structure is chosen to preserve compatibility with             !
  !    the current structure of PW/input.f90                               !
  !------------------------------------------------------------------------!
  SUBROUTINE read_xml_pw( attr, step )
    !
    !
    USE io_global, ONLY : ionode, ionode_id
    USE mp,        ONLY : mp_bcast
    USE iotk_module, ONLY : iotk_scan_begin, iotk_scan_end, iotk_scan_attr, iotk_attlenx
    USE iotk_unit_interf, ONLY : iotk_rewind
    !
    !
    IMPLICIT NONE
    !
    !
    INTEGER, INTENT(IN) :: step
    CHARACTER (len = *), INTENT(IN) :: attr
    !
    INTEGER :: ierr
    CHARACTER (len = iotk_attlenx) :: attr2
    CHARACTER (len = 30) :: field, card
    CHARACTER (len = 256) :: dummy
    LOGICAL :: found
    !
    !
    IF ( step==1 ) THEN  
       !
       ! ... reading of the attributes of the xml root node
       !
       CALL iotk_scan_attr( attr, 'calculation', dummy, found = found, ierr = ierr )
       IF ( .not. found ) CALL errore( 'read_xml_pw', 'attribute calculation of root node is &
            & compulsory', abs(ierr) )
       IF ( ierr /= 0 ) CALL errore( 'read_xml_pw', 'error reading calculation attribute of root node', 1 )
       calculation = trim( dummy )
       !
       !
       CALL iotk_scan_attr( attr, 'prefix', dummy, found = found, ierr = ierr )
       IF ( ierr /= 0 ) CALL errore( 'read_xml_pw', 'error reading prefix attribute of root node', abs(ierr) )
       IF ( found ) prefix = trim( dummy )
       !
       !
       CALL iotk_scan_attr( attr, 'title', dummy, found = found, ierr = ierr )
       IF ( ierr /= 0 ) CALL errore( 'read_xml_pw', 'error reading title attribute of root node', 1 )
       IF ( found ) title = trim( dummy )
       
       !
       ! ... fixing of some default values using the calculation variable
       !
       CALL fixval( 'PW' )
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
       !
       ! ... reading ATOMIC_LIST cards
       !
       CALL card_default( 'ATOMIC_LIST' )
       !
       IF ( ionode ) THEN
          !
          IF ( ( trim( calculation ) == 'neb' ) .or. ( trim( calculation ) == 'smd' ) ) THEN
             CALL card_xml_chain ( )
          ELSE
             CALL card_xml_atomic_list ( )
          END IF
          !
       END IF
       !
       CALL card_bcast( 'ATOMIC_LIST' )
       !
       !
       ! ... reading of all the FIELDS
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
       IF ( calculation /= 'nscf' .and. calculation /= 'bands'  ) THEN
          !
          IF ( restart_mode == 'from_scratch' ) THEN
             !
             startingwfc = 'atomic'
             startingpot = 'atomic'
             !
          ELSE
             !
             startingwfc = 'file'
             startingpot = 'file'
             !
          END IF
          !
       END IF
       
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
       CALL ee_bcast()
       CALL wannier_bcast()
       CALL wannier_ac_bcast()
       !
       !
       !
    ELSE
       ! ... second step : reading of the remaining cards
       !
       !
       ! ... reading of CONSTRAINTS cards
       !
       card = 'CONSTRAINTS'
       CALL card_default( card )
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
       CALL mp_bcast ( found, ionode_id )
       !
       IF ( found ) CALL card_bcast( 'CONSTRAINTS' )
       !
       !
       ! ... reading K_POINTS card
       !
       card = 'K_POINTS'
       CALL card_default( card )
       !
       IF ( ionode ) THEN
          !
          CALL iotk_scan_begin( xmlinputunit, trim( card ), attr = attr2, found = found, ierr = ierr )
          IF ( ierr /= 0 ) GO TO 9
          !
          IF ( found ) THEN
             !
             CALL card_xml_kpoints( attr2 )
             !
             CALL iotk_scan_end( xmlinputunit, trim( card ), ierr = ierr)
             IF ( ierr /= 0 ) GOTO 10
             !
          ELSE
             !
             GOTO 11
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
       card = 'OCCUPATIONS'
       CALL card_default( card )
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
       CALL mp_bcast ( found, ionode_id )
       !
       IF ( found ) CALL card_bcast( 'OCCUPATIONS' )
       !
       !
       ! ... reading CLIMBING_IMAGES card
       !
       card = 'CLIMBING_IMAGES'
       CALL card_default( card )
       !
       IF ( ionode ) THEN
          !
          CALL iotk_scan_begin( xmlinputunit, trim( card ), found = found, ierr = ierr )
          IF ( ierr /= 0 ) GO TO 9
          !
          IF ( found ) THEN
             !
             CALL card_xml_climbing_images()
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
       CALL mp_bcast ( found, ionode_id )
       !
       IF ( found ) CALL card_bcast( 'CLIMBING_IMAGES' )
       !
       !
    END IF
    !
    !
    RETURN
    !
9   CALL errore('read_xml_pw', 'error reading begin tag of '//card//' card', ABS( ierr ) )
10  CALL errore('read_xml_pw', 'error reading end tag of '//card//' card', ABS( ierr ) )
11  CALL errore('read_xml_pw', card//' card was not found', 1)
    !
    !
  END SUBROUTINE READ_XML_PW
  !
  !
  !
END MODULE read_xml_module
