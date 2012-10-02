!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM pwscf
  !----------------------------------------------------------------------------
  !
  ! ... Plane Wave Self-Consistent Field code 
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE cell_base,        ONLY : fix_volume, fix_area
  USE control_flags,    ONLY : conv_elec, gamma_only, lscf
  USE control_flags,    ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs
  USE force_mod,        ONLY : lforce, lstres, sigma
  USE environment,      ONLY : environment_start
  USE check_stop,       ONLY : check_stop_init, check_stop_now
  USE mp_global,        ONLY : mp_startup, mp_global_end, intra_image_comm
  USE mp_global,        ONLY : nimage, me_image, root_image, my_image_id
#if defined(__MS2)
  USE ms2,              ONLY : MS2_enabled,                 &
                               ms2_initialization,    &
                               set_positions, return_forces
#endif
  USE io_files,           ONLY : tmp_dir
  USE image_io_routines,  ONLY : io_image_start
  USE xml_io_base,        ONLY : create_directory, change_directory
  USE read_input,         ONLY : read_input_file
  !
  IMPLICIT NONE
  !
  !
  CHARACTER(len=256) :: dirname
  !
#ifdef __MPI
  !
  CALL mp_startup ( )
  ! reset IO nodes
  ! (do this to make each "image head node" an ionode)
  ! Has to be used ONLY to run nimage copies of pwscf
  !
  IF ( nimage > 1 ) CALL io_image_start( )
#endif
  CALL environment_start ( 'PWSCF' )
  !
  IF ( ionode ) WRITE( unit = stdout, FMT = 9010 ) &
         ntypx, npk, lmaxx
  !
  IF (ionode) CALL plugin_arguments()
  CALL plugin_arguments_bcast( ionode_id, intra_image_comm )
  !
  ! ... open, read, close input file
  !
  CALL read_input_file ('PW')
  !
  ! ... convert to internal variables
  !
  CALL iosys()
  !
  IF ( gamma_only ) WRITE( UNIT = stdout, &
     & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
  !
  IF( nimage > 1 ) THEN
     !
     ! ... When nimage are used, open a directory for each one
     ! ...It has to be done here in order not to disturb NEB like calculations
     !
     WRITE( dirname, FMT = '( I5.5 )' ) my_image_id
     tmp_dir = TRIM( tmp_dir )//TRIM( dirname )//'/'
     CALL create_directory( tmp_dir )
     CALL change_directory( tmp_dir )
     !
  END IF
  !
  ! call to void routine for user defined / plugin patches initializations
  !
  CALL plugin_initialization()
  !
  CALL check_stop_init()
  !
#if defined(__MS2)
  CALL ms2_initialization()
#endif
  !
  CALL setup ()
  !
#if defined(__MS2)
  CALL set_positions()
#endif
  !
  CALL init_run()
  !
  IF ( check_stop_now() ) CALL stop_run( .TRUE. )
  !
  main_loop: DO
     !
     ! ... electronic self-consistency
     !
     CALL electrons()
     !
     IF ( .NOT. conv_elec ) THEN
       CALL punch( 'all' )
       CALL stop_run( conv_elec )
     ENDIF
     !
     ! ... ionic section starts here
     !
     CALL start_clock( 'ions' )
     conv_ions = .TRUE.
     !
     ! ... recover from a previous run, if appropriate
     !
     IF ( restart .AND. lscf ) CALL restart_in_ions()
     !
     ! ... file in CASINO format written here if required
     !
     CALL pw2casino()
     !
     ! ... force calculation
     !
     IF ( lforce ) CALL forces()
     !
     ! ... stress calculation
     !
     IF ( lstres ) CALL stress ( sigma )
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        if (fix_volume) CALL impose_deviatoric_stress(sigma)
        !
        if (fix_area)  CALL  impose_deviatoric_stress_2d(sigma)
        !
        ! ... ionic step (for molecular dynamics or optimization)
        !
        CALL move_ions()
        !
        ! ... then we save restart information for the new configuration
        !
        IF ( istep < nstep .AND. .NOT. conv_ions ) THEN
           !
           CALL punch( 'config' )
           CALL save_in_ions()
           !
        END IF
        !
     END IF
     !
     CALL stop_clock( 'ions' )
     !
#if defined(__MS2)
     CALL return_forces()
#endif
     ! ... exit condition (ionic convergence) is checked here
     !
     IF ( conv_ions ) EXIT main_loop
     !
     ! ... terms of the hamiltonian depending upon nuclear positions
     ! ... are reinitialized here
     !
#if defined(__MS2)
     CALL set_positions()
#endif
     IF ( lmd .OR. lbfgs ) CALL hinit1()
     !
  END DO main_loop
  !
  ! ... save final data file
  !
  CALL punch('all')
  CALL stop_run( conv_ions )
  !
  !  END IF      
  !
  STOP
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !
END PROGRAM pwscf
