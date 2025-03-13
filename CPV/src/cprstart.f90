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
  !! Molecular Dynamics using Density-Functional Theory.  
  !! This is the main routine driver for Car-Parrinello simulations.  
  !! See the documentation coming with the Quantum ESPRESSO distribution
  !! for credits, references, appropriate citation of this code.
  !
  USE input,         ONLY : iosys_pseudo, iosys
  USE io_global,     ONLY : ionode, ionode_id
  USE environment,   ONLY : environment_start, print_cuda_info
  USE check_stop,    ONLY : check_stop_init
  USE mp_global,     ONLY : mp_startup
  USE mp_images,     ONLY : intra_image_comm
  USE read_input,    ONLY : read_input_file
  USE command_line_options, ONLY : input_file_
  !
  IMPLICIT NONE
  !
  ! ... program starts here
  !
  ! ... initialize MPI (parallel processing handling)
  !
  CALL mp_startup ( )
  !
  ! ... start the environment
  !
  CALL environment_start( 'CP' )
  CALL print_cuda_info() 
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
  CALL check_stop_init()
  !
  CALL cpr_loop( 1 )
  !
  CALL laxlib_end()
  CALL stop_cp_run()
  !
END PROGRAM main
