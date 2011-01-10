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
  USE io_global,        ONLY : ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE control_flags,    ONLY : conv_elec, conv_ions, lpath, gamma_only
  USE environment,      ONLY : environment_start, environment_end
  USE path_variables,   ONLY : conv_path
  USE check_stop,       ONLY : check_stop_init
  USE path_base,        ONLY : initialize_path, search_mep
  USE path_io_routines, ONLY : path_summary
  USE image_io_routines, ONLY : io_image_start
  USE mp_global,        ONLY : mp_startup, mp_bcast, mp_global_end
  !
  USE iotk_module,           ONLY : iotk_attlenx
  USE open_close_input_file_interf, ONLY : open_input_file, close_input_file
  USE read_xml_module,       ONLY : read_xml
  USE read_cards_module,     ONLY : read_cards
  USE read_namelists_module, ONLY : read_namelists
  USE path_read_namelists_module, ONLY : path_read_namelist
  USE path_read_cards_module, ONLY : path_read_cards
  !
  USE path_io_units_module, ONLY : stdinpath, set_input_unit
  !
  USE path_input_parameters_module, ONLY : nstep_path
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
  !
  call set_input_unit()
  !
  open(unit=stdinpath,file="neb.dat",status="old")
  CALL path_read_namelist(stdinpath)
  !
  IF ( xmlinput ) THEN
     CALL read_xml ('PW', attr = attr )
  ELSE
     CALL read_namelists( prog='PW', unit=5 )

     CALL set_engine_input_defaults()

     CALL read_cards( prog='PW', unit=5 )
  ENDIF
  !
  CALL path_read_cards(stdinpath)
  !
  !
  close(stdinpath)
  !
  CALL set_engine_io_units()
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
  !
  STOP
  !
9010 FORMAT( /,5X,'Current dimensions of program SM are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !
END PROGRAM sm
