!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE run_pwscf ( conv ) 
  !----------------------------------------------------------------------------
  !
  ! ... Run an instance of the Plane Wave Self-Consistent Field code 
  ! ... MPI initialization and input data reading is performed in the 
  ! ... calling code - returns in "conv" whether successfully completed
  ! ... Will be eventualy merged with NEB
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE cell_base,        ONLY : fix_volume, fix_area
  USE control_flags,    ONLY : conv_elec, gamma_only, lscf
  USE control_flags,    ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs
  USE force_mod,        ONLY : lforce, lstres, sigma
  USE check_stop,       ONLY : check_stop_init, check_stop_now
  USE mp_images,        ONLY : intra_image_comm
#if defined(__MS2)
  USE ms2,              ONLY : MS2_enabled,                 &
                               ms2_initialization,    &
                               set_positions, return_forces
#endif
  !
  IMPLICIT NONE
  LOGICAL, INTENT(OUT) :: conv
  !
  !
  IF ( ionode ) WRITE( unit = stdout, FMT = 9010 ) ntypx, npk, lmaxx
  !
  IF (ionode) CALL plugin_arguments()
  CALL plugin_arguments_bcast( ionode_id, intra_image_comm )
  !
  ! ... convert to internal variables
  !
  CALL iosys()
  !
  IF ( gamma_only ) WRITE( UNIT = stdout, &
     & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
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
  IF ( check_stop_now() ) THEN
     CALL punch( 'all' )
     conv = .TRUE.
     RETURN
  ENDIF
  !
  main_loop: DO
     !
     ! ... electronic self-consistency
     !
     CALL electrons()
     !
     IF ( .NOT. conv_elec ) THEN
       CALL punch( 'all' )
       conv = .FALSE.
       RETURN
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
     IF ( lmd ) CALL pw2casino()
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
  IF ( .not. lmd) CALL pw2casino()
  CALL punch('all')
  !
  conv = conv_ions
  RETURN
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !
END SUBROUTINE run_pwscf
