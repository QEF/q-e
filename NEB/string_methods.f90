!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM sm
  !----------------------------------------------------------------------------
  !
  ! ... Plane Wave Self-Consistent Field code 
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE control_flags,    ONLY : conv_elec, conv_ions, lpath, gamma_only
  USE environment,      ONLY : environment_start, environment_end
  USE path_variables,   ONLY : conv_path
  USE check_stop,       ONLY : check_stop_init
  USE path_base,        ONLY : initialize_path, search_mep
  USE path_io_routines, ONLY : path_summary
  USE image_io_routines, ONLY : io_image_start
  USE mp_global,        ONLY : mp_startup, mp_bcast, mp_global_end
  USE read_namelists_module, ONLY : read_namelists
  !
  USE iotk_module,           ONLY : iotk_attlenx
  USE open_close_input_file_interf, ONLY : open_input_file, close_input_file
  USE read_xml_module,       ONLY : read_xml
  !
  IMPLICIT NONE
  !
  !
  LOGICAL :: xmlinput = .false.
  CHARACTER (len=iotk_attlenx) :: attr
  !
  !
#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'SM' )
  !
  IF ( ionode ) THEN
     !
#if defined (EXX)
     WRITE( UNIT = stdout, &
         & FMT = '(/,5X,"EXPERIMENTAL VERSION WITH EXACT EXCHANGE")' )
#endif
     WRITE( unit = stdout, FMT = 9010 ) &
         ntypx, npk, lmaxx
     !
  END IF   
  !
  ! INPUT RELATED
  !
  ! ... open input file
  !
  IF( ionode ) CALL open_input_file(xmlinput,attr)
  !
  ! bcast of xmlinput and attr needs to be done
  ! because is only the open statement inside
  ! read_cards and read_namelist (in Modules) that has
  ! if(ionode) !!! in future call read_cards_pw, call read_namelis
  ! call read_xml should be done only by ionode. bcast is already done
  ! inside read_cards and read_namelist.
  !
  call mp_bcast(xmlinput,ionode_id)
  call mp_bcast(attr,ionode_id)
  !
  IF( xmlinput ) THEN
    CALL read_xml( 'PW', attr )
  ELSE 
    CALL read_namelists('SM')
  ENDIF
  !
  CALL set_defaults()
  !
  CALL iosys(xmlinput,attr)
  !
  !
  CALL ioneb(xmlinput,attr)
  ! ... close_input_file(xmlinput)
  !
  IF( ionode ) CALL close_input_file(xmlinput)
  !
  ! END INPUT RELATED
  !
  !
  CALL check_stop_init()
  !
  CALL io_image_start()
  !
  CALL initialize_path()
  !
  CALL path_summary()
  !
  CALL search_mep()
  !
  CALL stop_run_path( conv_path )
  !
!  CALL environment_end( 'SM' )
  !
!  CALL mp_global_end()
  !
  STOP
  !
9010 FORMAT( /,5X,'Current dimensions of program SM are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !
END PROGRAM sm
