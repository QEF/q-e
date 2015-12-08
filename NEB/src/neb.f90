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
  USE io_global,         ONLY : meta_ionode, meta_ionode_id
  USE environment,       ONLY : environment_start, environment_end
  USE check_stop,        ONLY : check_stop_init
  USE mp,                ONLY : mp_bcast
  USE mp_global,         ONLY : mp_startup
  USE mp_world,          ONLY : world_comm, mpime, root
  USE read_input,        ONLY : read_input_file
  USE command_line_options,  ONLY : input_file_
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
  USE path_io_units_module,  ONLY : iunpath
  !
  IMPLICIT NONE
  !
  CHARACTER(len=256) :: engine_prefix, parsing_file_name
  INTEGER :: unit_tmp, i, iimage
  INTEGER, EXTERNAL :: find_free_unit, input_images_getarg
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  !
  CALL mp_startup ( start_images=.true. )
  CALL environment_start ( 'NEB' )
  !
  ! INPUT RELATED
  !
  engine_prefix = "pw_"
  !
  IF (meta_ionode) CALL plugin_arguments()
  CALL plugin_arguments_bcast( meta_ionode_id, world_comm )
  !
  ! ... open input file
  !
  IF ( input_file_ /= ' ') THEN
     WRITE(iunpath,'(/,5X,"parsing_file_name: ",A)') trim(input_file_)
     CALL path_gen_inputs ( trim(input_file_), engine_prefix, &
                            input_images, root, world_comm )
  ELSE
     WRITE(iunpath,'(/,5X,"No input file found, assuming nothing to parse",/,&
    &               5X,"Searching argument -input_images or --input_images")')
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
  do i=1,input_images
    !
    IF ( i > 1 ) CALL clean_pw(.true.)
    parsing_file_name = trim(engine_prefix)//trim(int_to_char(i))//".in"
    !
    CALL read_input_file ( 'PW', parsing_file_name )
    CALL iosys()
    !
    IF ( i == 1 ) THEN
      CALL engine_to_path_nat()
      CALL engine_to_path_alat()
      CALL allocate_path_input_ions(input_images)
    END IF
    CALL engine_to_path_pos(i)
    IF ( i == 1 ) CALL engine_to_path_fix_atom_pos()
    !
  enddo
  !
  CALL path_to_engine_fix_atom_pos()
  !
  CALL ioneb()
  CALL set_engine_output()
  !
  ! END INPUT RELATED
  !
  CALL check_stop_init()
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
