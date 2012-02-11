!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM neb
  !----------------------------------------------------------------------------
  !
  ! ... Nudged Elastic Band / Strings Method algorithm
  !
  USE io_global,        ONLY : meta_ionode_id, xmlinputunit
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE control_flags,    ONLY : conv_elec, conv_ions, lpath, gamma_only
  USE environment,      ONLY : environment_start, environment_end
  USE path_variables,   ONLY : conv_path
  USE check_stop,       ONLY : check_stop_init
  USE path_base,        ONLY : initialize_path, search_mep
  USE path_io_routines, ONLY : path_summary
  USE image_io_routines, ONLY : io_image_start
  USE mp_global,         ONLY : mp_bcast, mp_rank
  !
  USE mp_image_global_module, ONLY : mp_image_startup, world_comm
  USE mp_image_global_module, ONLY : me_image, nimage
  USE mp_global,             ONLY : mp_start
  USE iotk_module,           ONLY : iotk_attlenx
  USE open_close_input_file, ONLY : open_input_file, close_input_file
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
  USE iotk_module,   ONLY : iotk_open_read, iotk_close_read,iotk_attlenx
  !
  IMPLICIT NONE
  !
  !
  CHARACTER (len=iotk_attlenx) :: attr
  !
  CHARACTER(len=256) :: engine_prefix
  !
  INTEGER :: unit_tmp
  !
  INTEGER :: i, iimage
  CHARACTER(len=10) :: a_tmp
  !
  INTEGER :: mpime = 0, nproc = 1, neb_comm = 0
  INTEGER :: root = 0 
  !
  CHARACTER(len=256) :: parsing_file_name
  LOGICAL :: lfound_parsing_file, lfound_input_images
  !
  LOGICAL :: lxml
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  unit_tmp = 45
  xmlinputunit = 45
  !
#ifdef __MPI
  CALL mp_start(nproc,mpime,neb_comm)
  CALL mp_image_startup (root,neb_comm)
  CALL engine_mp_start()
#endif
  CALL environment_start ( 'NEB' )
  !
  !
  ! INPUT RELATED
  !
  ! ... open input file
  !
  if(mpime==root) CALL input_file_name_getarg(parsing_file_name,lfound_parsing_file)
  !
  engine_prefix = "pw_"
  !
  CALL mp_bcast(parsing_file_name,root,neb_comm)
  CALL mp_bcast(lfound_parsing_file,root,neb_comm)
  !
  !
  if(lfound_parsing_file) then
    write(0,*) ""
    write(0,*) "parsing_file_name: ", trim(parsing_file_name)
  call path_gen_inputs(trim(parsing_file_name),engine_prefix,input_images,root,neb_comm)
  !
  else
  !
    write(0,*) ""
    write(0,*) "NO input file found, assuming nothing to parse."
    write(0,*) "Searching argument -input_images or --input_images"
    if(mpime==root) CALL input_images_getarg(input_images,lfound_input_images)
    CALL mp_bcast(input_images,root,neb_comm)
    CALL mp_bcast(lfound_input_images,root,neb_comm)
  !
    IF(.not.lfound_input_images) CALL errore('string_methods', &
      'Neither a file to parse nor input files for each image found',1)
  !
  endif
  !
! mpi for engine
!  call path_to_engine_mp()
  !  
  !
  call set_input_unit()
  !
  open(unit=stdinpath,file="neb.dat",status="old")
  CALL path_read_namelist(stdinpath)
  CALL path_read_cards(stdinpath)
  close(stdinpath)
  !
  !
  OPEN(unit_tmp, file=trim(engine_prefix)//"1.in")
  CALL test_input_xml(unit_tmp,lxml)
  CLOSE(unit_tmp)
  if(.not.lxml) then
    OPEN(unit_tmp, file=trim(engine_prefix)//"1.in")
    CALL read_namelists( prog='PW', unit=unit_tmp )
    CALL read_cards( prog='PW', unit=unit_tmp )
    CLOSE(unit_tmp)
  else
    CALL iotk_open_read( unit_tmp, trim(engine_prefix)//"1.in", &
           attr = attr, qe_syntax = .true.)
    CALL read_xml('PW', attr = attr )
    CALL iotk_close_read(unit_tmp)
  endif
  CALL iosys()
  CALL engine_to_path_nat()
  CALL engine_to_path_alat()
  CALL allocate_path_input_ions(input_images)
  CALL engine_to_path_pos(1)
  CALL engine_to_path_fix_atom_pos()
  do i=2,input_images
    CALL set_engine_input_defaults()
    CALL clean_pw(.true.)
!    if(i>=1.and.i<10) then
!    write(a_tmp,'(i1)') i
!    elseif(i>10.and.i<100) then
!    write(a_tmp,'(i2)') i
!    elseif(i>100.and.i<1000) then
!    write(a_tmp,'(i3)')
!    endif

    a_tmp=trim(int_to_char(i))

    OPEN(unit_tmp,file=trim(engine_prefix)//trim(a_tmp)//".in") 
    CALL test_input_xml(unit_tmp,lxml)
    CLOSE(unit_tmp)
    if(.not.lxml) then
      OPEN(unit_tmp,file=trim(engine_prefix)//trim(a_tmp)//".in")
      CALL read_namelists( prog='PW', unit=unit_tmp )
      CALL read_cards( prog='PW', unit=unit_tmp )
      CLOSE(unit_tmp)
    else
      CALL iotk_open_read( unit_tmp, trim(engine_prefix)//trim(a_tmp)//".in", &
           attr = attr, qe_syntax = .true. )
      CALL read_xml('PW', attr = attr )
      CALL iotk_close_read(unit_tmp)
    endif
    CALL iosys()
    CALL engine_to_path_pos(i)
  enddo
  !
  !
  CALL path_to_engine_fix_atom_pos()
  !
  CALL ioneb()
  !
  !
  CALL set_engine_io_units()
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
END PROGRAM neb
