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
  USE io_global,        ONLY : meta_ionode_id
  USE io_files,         ONLY : find_free_unit
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE control_flags,    ONLY : conv_elec, conv_ions, lpath, gamma_only
  USE environment,      ONLY : environment_start, environment_end
  USE path_variables,   ONLY : conv_path
  USE check_stop,       ONLY : check_stop_init
  USE path_base,        ONLY : initialize_path, search_mep
  USE path_io_routines, ONLY : path_summary
  USE image_io_routines, ONLY : io_image_start
  USE mp_global,        ONLY : mp_startup, mp_bcast, mp_global_end, world_comm
  !
  USE iotk_module,           ONLY : iotk_attlenx
  USE open_close_input_file_interf, ONLY : open_input_file, close_input_file
  USE read_xml_module,       ONLY : read_xml
  USE read_cards_module,     ONLY : read_cards
  USE read_namelists_module, ONLY : read_namelists
  USE path_read_namelists_module, ONLY : path_read_namelist
  USE path_read_cards_module,            ONLY : path_read_cards
  !
  USE path_io_units_module, ONLY : stdinpath, set_input_unit
  !
  USE path_input_parameters_module, ONLY : nstep_path
  !
  USE path_input_parameters_module, ONLY : input_images
  !
  USE path_input_parameters_module, ONLY : allocate_path_input_ions
  !
  IMPLICIT NONE
  !
  !
  LOGICAL :: xmlinput = .false.
  CHARACTER (len=iotk_attlenx) :: attr
  !
  CHARACTER(len=256) :: engine_prefix
  !
  INTEGER :: unit_tmp
  !
  INTEGER :: i
  CHARACTER(len=10) :: a_tmp
  !
  !
  unit_tmp = 45
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
  !
  engine_prefix = "pw_"
  call path_gen_inputs("myinput.in",engine_prefix,input_images,meta_ionode_id,world_comm)
  !
  call set_input_unit()
  !
  open(unit=stdinpath,file="neb.dat",status="old")
  CALL path_read_namelist(stdinpath)
  CALL path_read_cards(stdinpath)
  close(stdinpath)
  !
  !
  CALL set_engine_io_units()

!  unit_tmp = find_free_unit()
  OPEN(unit_tmp, file=trim(engine_prefix)//"1.in")
  CALL read_namelists( prog='PW', unit=unit_tmp )
  CALL read_cards( prog='PW', unit=unit_tmp )
  CALL iosys()
  CALL engine_to_path_nat()
  CALL engine_to_path_alat()
  CALL allocate_path_input_ions(input_images)
  CALL engine_to_path_pos(1)
  CALL engine_to_path_fix_atom_pos()
  CLOSE(unit_tmp)

  do i=2,input_images
    CALL set_engine_input_defaults()
    CALL clean_pw(.true.)
    if(i>=1.and.i<10) then
    write(a_tmp,'(i1)') i
    elseif(i>10.and.i<100) then
    write(a_tmp,'(i2)') i
    elseif(i>100.and.i<1000) then
    write(a_tmp,'(i3)')
    endif

!    unit_tmp = find_free_unit()
    OPEN(unit_tmp,file=trim(engine_prefix)//trim(a_tmp)//".in") 
    CALL read_namelists( prog='PW', unit=unit_tmp )
    CALL read_cards( prog='PW', unit=unit_tmp )
    CALL iosys()
    CALL engine_to_path_pos(i)
    CLOSE(unit_tmp)
  enddo
  !
  !
  CALL path_to_engine_fix_atom_pos()
  !
  CALL ioneb(xmlinput,attr)
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
