!
! Copyright (C) 2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Aadhityan A, Lorenzo Paulatto, Michele Casula, Tommaso Morresi
!
!----------------------------------------------------------------------------
PROGRAM pioud
  !----------------------------------------------------------------------------
  !
  ! ... PATH INTEGRAL MOLECULAR DYNAMICS CODE !!!
  !
  USE io_global,                    ONLY : meta_ionode, meta_ionode_id
  USE environment,                  ONLY : environment_start, environment_end
  USE check_stop,                   ONLY : check_stop_init
  USE mp,                           ONLY : mp_bcast
  USE mp_global,                    ONLY : mp_startup
  USE mp_world,                     ONLY : world_comm, mpime, root
  USE mp_pools,                     ONLY : intra_pool_comm
  USE mp_bands,                     ONLY : intra_bgrp_comm, inter_bgrp_comm
  USE read_input,                   ONLY : read_input_file
  USE command_line_options,         ONLY : input_file_, ndiag_
  USE trpmd_base,                   ONLY : initialize_polymer, explore_phasespace
  USE ring_input_parameters_module, ONLY : input_images, allocate_path_input_ions, &
                                           deallocate_path_input_ions
  USE ring_io_units_module,         ONLY : iunpath

  IMPLICIT NONE
  !
  !include 'laxlib.fh'
  !
  CHARACTER(len=256) :: engine_prefix, parsing_file_name
  INTEGER :: unit_tmp, i, iimage
  INTEGER, EXTERNAL :: find_free_unit, input_images_getarg
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  
  CALL mp_startup ( start_images=.true.)
  !
  CALL environment_start ( 'PIOUD' )   !!! <----my mod.
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
     
     CALL pimd_gen_inputs( trim(input_file_), engine_prefix, &                 !!! <----my mod.
                            root, world_comm  )                                !!! <----my mod.
     CALL pw_gen_inputs( trim(input_file_), engine_prefix, & 
                            input_images, root, world_comm ) !!! Contains NEB. But can't get ridoff completly since it require to generate pw_1.in 
  !Remove Else or just put no input images.
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
  

  
  IF ( meta_ionode) THEN                  !!! <----my mod.
    !  unit_tmp = find_free_unit ()                  !!! <----my mod.
     open(newunit=unit_tmp,file="pimd.dat",status="old")      !!! <----my mod.
     CALL pimd_read_input(unit_tmp)                !!! <----my mod.
     close(unit=unit_tmp)                 !!! <----my mod.
  END IF                                  !!! <----my mod.
  

  
  CALL pimd_mp_bcast   !!! <----my mod.
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
      CALL engine_to_path_nat_alat()
      CALL allocate_path_input_ions(input_images)
    END IF
    CALL engine_to_path_pos(i)
    !
  enddo
  
  CALL pimd_get_amas_and_nat

  CALL verify_pioud_tmpdir()
  CALL set_engine_output()

  CALL check_stop_init() 
  CALL initialize_polymer()
  CALL deallocate_path_input_ions()

  CALL explore_phasespace()

  CALL laxlib_end()
  
  CALL pimd_deallocation    !!! <----my mod.
  if ( meta_ionode) call pimd_close_files   !!! <----my mod.
  
  CALL stop_run_path( .true. )     !!! <----my mod.
  
  STOP
  
END PROGRAM pioud
