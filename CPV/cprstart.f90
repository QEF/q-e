!
! Copyright (C) 2002-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!==============================================================================
!***  Molecular Dynamics using Density-Functional Theory                   ****
!***  this is the main routine driver for Car-Parrinello simulations       ****
!******************************************************************************
!***  See the documentation coming with the Quantum-Espresso distribution  ****
!***  for credits, references, appropriate citation of this code           ****
!******************************************************************************
!
!----------------------------------------------------------------------------
PROGRAM main
  !----------------------------------------------------------------------------
  !
  USE input,         ONLY : read_input_file, iosys_pseudo, iosys
  USE io_global,     ONLY : io_global_start, io_global_getmeta
  USE mp_global,     ONLY : mp_global_start, init_pool
  USE mp,            ONLY : mp_end, mp_start, mp_env, mp_bcast
  USE control_flags, ONLY : lneb, lsmd, lmetadyn, program_name
  USE control_flags, ONLY : use_task_groups
  USE environment,   ONLY : environment_start
  USE check_stop,    ONLY : check_stop_init
  !
  IMPLICIT NONE
  !
  INTEGER            :: mpime, nproc, world, meta_ionode_id
  INTEGER            :: nimage, ntask_groups
  LOGICAL            :: meta_ionode
  INTEGER, PARAMETER :: root = 0
  !
  ! ... program starts here
  !
  program_name = 'CP90'
  !
  ! ... initialize MPI (parallel processing handling)
  !
  CALL mp_start()
  !
  ! ... get from communication sub-sistem basic parameters
  ! ... to handle processors
  !
  CALL mp_env( nproc, mpime, world )
  !
  ! ... now initialize module holding processors and groups
  ! ... variables
  !
  CALL mp_global_start( root, mpime, world, nproc )
  !
  ! ... mpime = processor number, starting from 0
  ! ... nproc = number of processors
  ! ... world = group index of all processors
  ! ... root  = index of the root processor
  !
  ! ... initialize input output
  !
  CALL io_global_start( mpime, root )
  !
  ! ... get the "meta" io node
  !
  CALL io_global_getmeta( meta_ionode, meta_ionode_id )
  !
  IF ( meta_ionode ) THEN
     !
     ! ... check for command line arguments
     !
     CALL get_arg_nimage( nimage )
     !
     nimage = MAX( nimage, 1 )
     nimage = MIN( nimage, nproc )
     !
     CALL get_arg_ntg( ntask_groups )
     !
  END IF
  !
  CALL mp_bcast( nimage,       meta_ionode_id, world )
  CALL mp_bcast( ntask_groups, meta_ionode_id, world )
  !
  IF( ntask_groups > 1 ) THEN
     use_task_groups = .TRUE.
  END IF
  !
  ! ... start the environment
  !
  CALL environment_start( )
  !
  ! ... readin the input file
  !
  CALL read_input_file()
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
  CALL check_stop_init()
  !
  ! ... here reorganize processors in groups
  !
  CALL init_pool( nimage, ntask_groups )
  !
  IF ( lneb ) THEN
     !
     CALL neb_loop( )
     !
  ELSE IF ( lsmd ) THEN
     !
     CALL errore ( 'cpr_main', 'SMD no longer implemented', 1)
     !
  ELSE
     !
     CALL cpr_loop( 1 )
     !
  END IF
  !
  CALL stop_run( .TRUE. )
  !
  STOP
  !
END PROGRAM main
