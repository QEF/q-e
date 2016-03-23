!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!==============================================================================
!***  Molecular Dynamics using Density-Functional Theory                   ****
!***  this is the main routine driver for Car-Parrinello simulations       ****
!******************************************************************************
!***  See the documentation coming with the Quantum ESPRESSO distribution  ****
!***  for credits, references, appropriate citation of this code           ****
!******************************************************************************
!
!----------------------------------------------------------------------------
PROGRAM main
  !----------------------------------------------------------------------------
  !
  USE input,         ONLY : iosys_pseudo, iosys
  USE read_input,    ONLY : read_input_file
  USE mp_global,     ONLY : mp_startup
  USE io_global,     ONLY : ionode, ionode_id
  USE environment,   ONLY : environment_start
  USE check_stop,    ONLY : check_stop_init
  USE mp_images,     ONLY : intra_image_comm
  USE command_line_options, ONLY : input_file_
  !
  IMPLICIT NONE
  !
  ! ... program starts here
  !
  ! ... initialize MPI (parallel processing handling)
  !
  CALL mp_startup ( diag_in_band_group = .true. )
  !
  ! ... start the environment
  !
  CALL environment_start( 'CP' )
  !
  ! reading plugin arguments
  !
  IF(ionode) CALL plugin_arguments()
  CALL plugin_arguments_bcast(ionode_id,intra_image_comm)
  !
  ! ... open, read, close the input file
  !
  CALL read_input_file( 'CP', input_file_ )
  !
  ! ... read in pseudopotentials files and then
  ! ... copy pseudopotential parameters into internal variables
  !
  CALL iosys_pseudo()
  !
  ! ... copy-in input parameters from input_parameter module
  !
  CALL iosys()
  !
  ! call to void routine for user define / plugin patches initializations
  ! temporary moved to init_run
!  CALL plugin_initialization()
  !
  !
  CALL check_stop_init()
  !
  CALL cpr_loop( 1 )
  !
  CALL stop_run()
  CALL do_stop( .TRUE. )
  !
  STOP
  !
END PROGRAM main
