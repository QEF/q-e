!
! Copyright (C) 2011-2013 Quantum ESPRESSO group
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
  USE io_global,         ONLY : meta_ionode_id
  USE environment,       ONLY : environment_start, environment_end
  USE check_stop,        ONLY : check_stop_init
  USE image_io_routines, ONLY : io_image_start
  USE mp,                ONLY : mp_bcast
  USE mp_global,         ONLY : mp_startup, nimage, world_comm, mpime, root
  USE iotk_module,       ONLY : iotk_open_read, iotk_close_read, iotk_attlenx
  USE read_xml_module,       ONLY : read_xml
  USE read_cards_module,     ONLY : read_cards
  USE read_namelists_module, ONLY : read_namelists
  !
  USE path_variables,    ONLY : conv_path
  USE path_base,         ONLY : initialize_path, search_mep
  USE path_io_routines,  ONLY : path_summary
  USE path_read_namelists_module, ONLY : path_read_namelist
  USE path_read_cards_module,     ONLY : path_read_cards
  !
  USE path_input_parameters_module, ONLY : nstep_path, input_images, &
                                           allocate_path_input_ions, &
                                           deallocate_path_input_ions
  !
  IMPLICIT NONE
  !
  CHARACTER (len=iotk_attlenx) :: attr
  CHARACTER(len=256) :: engine_prefix
  !
  INTEGER :: unit_tmp
  INTEGER :: i, iimage
  CHARACTER(len=10) :: a_tmp
  !
  CHARACTER(len=256) :: parsing_file_name
  LOGICAL :: lfound_parsing_file, lxml
  !
  INTEGER, EXTERNAL :: find_free_unit, input_images_getarg
  LOGICAL, EXTERNAL :: test_input_xml, input_file_name_getarg
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  !
#ifdef __MPI
  CALL mp_startup ( start_images=.true. )
  IF ( nimage > 1 ) CALL io_image_start( )
#endif
  CALL environment_start ( 'NEB' )
  !
  ! INPUT RELATED
  !
  engine_prefix = "pw_"
  !
  ! ... open input file
  !
  IF ( mpime == root) &
     lfound_parsing_file = input_file_name_getarg (parsing_file_name)
  CALL mp_bcast(parsing_file_name,root,world_comm)
  CALL mp_bcast(lfound_parsing_file,root,world_comm)
  !
  IF (lfound_parsing_file) then
     WRITE(0,*) ""
     WRITE(0,*) "parsing_file_name: ", trim(parsing_file_name)
     CALL path_gen_inputs ( trim(parsing_file_name), engine_prefix, &
                            input_images, root, world_comm )
  ELSE
     WRITE(0,*) ""
     WRITE(0,*) "NO input file found, assuming nothing to parse."
     WRITE(0,*) "Searching argument -input_images or --input_images"
     IF ( mpime == root )  input_images = input_images_getarg ( )
     CALL mp_bcast(input_images,root, world_comm)
     !
     IF (input_images == 0) CALL errore('string_methods', &
        'Neither a file to parse nor input files for each image found',1)
     !
  ENDIF
  !
  unit_tmp = find_free_unit () 
  open(unit=unit_tmp,file="neb.dat",status="old")
  CALL path_read_namelist(unit_tmp)
  CALL path_read_cards(unit_tmp)
  close(unit=unit_tmp)
  !
  OPEN(unit_tmp, file=trim(engine_prefix)//"1.in")
  lxml = test_input_xml(unit_tmp)
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

    a_tmp=trim(int_to_char(i))
    OPEN(unit_tmp,file=trim(engine_prefix)//trim(a_tmp)//".in") 
    lxml = test_input_xml(unit_tmp)
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
    !
    CALL iosys()
    !
    CALL engine_to_path_pos(i)
    !
  enddo
  !
  CALL path_to_engine_fix_atom_pos()
  !
  CALL ioneb()
  CALL set_engine_io_units()
  !
  ! END INPUT RELATED
  !
  CALL check_stop_init()
  CALL io_image_start()
  CALL initialize_path()
  CALL deallocate_path_input_ions()
  CALL path_summary()
  !
  CALL search_mep()
  !
  CALL stop_run_path( conv_path )
  !
  STOP
  !
END PROGRAM neb
