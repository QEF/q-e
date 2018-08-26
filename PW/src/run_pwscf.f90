!
! Copyright (C) 2013-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE run_pwscf ( exit_status ) 
  !----------------------------------------------------------------------------
  !
  !! author: Paolo Giannozzi
  !! license: GNU 
  !! summary: Run an instance of the Plane Wave Self-Consistent Field code
  !!
  !! Run an instance of the Plane Wave Self-Consistent Field code 
  !! MPI initialization and input data reading is performed in the 
  !! calling code - returns in exit_status the exit code for pw.x, 
  !! returned in the shell. Values are:
  !! * 0: completed successfully
  !! * 1: an error has occurred (value returned by the errore() routine)
  !! * 2-127: convergence error
  !!   * 2: scf convergence error
  !!   * 3: ion convergence error
  !! * 128-255: code exited due to specific trigger
  !!   * 255: exit due to user request, or signal trapped,
  !!          or time > max_seconds
  !!     (note: in the future, check_stop_now could also return a value
  !!     to specify the reason of exiting, and the value could be used
  !!     to return a different value for different reasons)
  !! @Note
  !! 10/01/17 Samuel Ponce: Add Ford documentation
  !! @endnote
  !!
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE parameters,       ONLY : ntypx, npk, lmaxx
  USE cell_base,        ONLY : fix_volume, fix_area
  USE control_flags,    ONLY : conv_elec, gamma_only, ethr, lscf, twfcollect
  USE control_flags,    ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs
  USE command_line_options, ONLY : command_line
  USE force_mod,        ONLY : lforce, lstres, sigma, force
  USE check_stop,       ONLY : check_stop_init, check_stop_now
  USE mp_images,        ONLY : intra_image_comm
  USE extrapolation,    ONLY : update_file, update_pot
  USE scf,              ONLY : rho
  USE lsda_mod,         ONLY : nspin
  USE fft_base,         ONLY : dfftp
  USE qmmm,             ONLY : qmmm_initialization, qmmm_shutdown, &
                               qmmm_update_positions, qmmm_update_forces
  USE qexsd_module,     ONLY : qexsd_set_status
  !
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: exit_status
  !! Gives the exit status at the end
  LOGICAL, external :: matches
  !! checks if first string is contained in the second
  INTEGER :: idone 
  !! counter of electronic + ionic steps done in this run
  INTEGER :: ions_status = 3
  !!    ions_status =  3  not yet converged
  !!    ions_status =  2  converged, restart with nonzero magnetization
  !!    ions_status =  1  converged, final step with current cell needed
  !!    ions_status =  0  converged, exiting
  !
  exit_status = 0
  IF ( ionode ) WRITE( unit = stdout, FMT = 9010 ) ntypx, npk, lmaxx
  !
  IF (ionode) CALL plugin_arguments()
  CALL plugin_arguments_bcast( ionode_id, intra_image_comm )
  !
  ! ... needs to come before iosys() so some input flags can be
  !     overridden without needing to write PWscf specific code.
  ! 
  CALL qmmm_initialization()
  !
  ! ... convert to internal variables
  !
  CALL iosys()
  !
  ! ... If executable names is "dist.x", compute atomic distances, angles,
  ! ... nearest neighbors, write them to file "dist.out", exit
  !
  IF ( matches('dist.x',command_line) ) THEN
     IF (ionode) CALL run_dist ( exit_status )
     RETURN
  END IF
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
  CALL setup ()
  !
  CALL qmmm_update_positions()
  !
  ! ... dry run: code will stop here if called with exit file present
  ! ... useful for a quick and automated way to check input data
  !
  IF ( check_stop_now() ) THEN
     CALL pre_init()
     CALL data_structure( gamma_only )
     CALL summary()
     CALL memory_report()
     CALL qexsd_set_status(255)
     CALL punch( 'config' )
     exit_status = 255
     RETURN
  ENDIF
  !
  CALL init_run()
  !
  IF ( check_stop_now() ) THEN
     CALL qexsd_set_status(255)
     CALL punch( 'config' )
     exit_status = 255
     RETURN
  ENDIF
  !
  main_loop: DO idone = 1, nstep
     !
     ! ... electronic self-consistency or band structure calculation
     !
     IF ( .NOT. lscf) THEN
        CALL non_scf ()
     ELSE
        CALL electrons()
     END IF
     !
     ! ... code stopped by user or not converged
     !
     IF ( check_stop_now() .OR. .NOT. conv_elec ) THEN
        IF ( check_stop_now() ) exit_status = 255
        IF ( .NOT. conv_elec )  exit_status =  2
        CALL qexsd_set_status(exit_status)
        ! workaround for the case of a single k-point
        twfcollect = .FALSE.
        CALL punch( 'config' )
        RETURN
     ENDIF
     !
     ! ... ionic section starts here
     !
     CALL start_clock( 'ions' ); !write(*,*)' start ions' ; FLUSH(6)
     conv_ions = .TRUE.
     !
     ! ... recover from a previous run, if appropriate
     !
     !IF ( restart .AND. lscf ) CALL restart_in_ions()
     !
     ! ... file in CASINO format written here if required
     !
     IF ( lmd ) THEN
        CALL pw2casino( istep )
     ELSE
        CALL pw2casino( 0 )
     END IF
     !
     ! ... force calculation
     !
     IF ( lforce ) CALL forces()
     !
     ! ... stress calculation
     !
     IF ( lstres ) CALL stress ( sigma )
     !
     ! ... send out forces to MM code in QM/MM run
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        if (fix_volume) CALL impose_deviatoric_stress(sigma)
        if (fix_area)  CALL  impose_deviatoric_stress_2d(sigma)
        !
        ! ... save data needed for potential and wavefunction extrapolation
        !
        CALL update_file ( )
        !
        ! ... ionic step (for molecular dynamics or optimization)
        !
        CALL move_ions ( idone, ions_status )
        conv_ions = ( ions_status == 0 )
        !
        ! ... then we save restart information for the new configuration
        !
        IF ( idone <= nstep .AND. .NOT. conv_ions ) THEN 
            CALL qexsd_set_status(255)
            CALL punch( 'config' )
        END IF
        !
     END IF
     !
     CALL stop_clock( 'ions' ); !write(*,*)' stop ions' ; FLUSH(6)
     !
     CALL qmmm_update_forces( force, rho%of_r, nspin, dfftp)
     !
     ! ... exit condition (ionic convergence) is checked here
     !
     IF ( lmd .OR. lbfgs ) CALL add_qexsd_step(idone)
     IF ( conv_ions ) EXIT main_loop
     !
     ! ... receive new positions from MM code in QM/MM run
     !
     CALL qmmm_update_positions()
     !
     ! ... terms of the hamiltonian depending upon nuclear positions
     ! ... are reinitialized here
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        IF ( ions_status == 1 ) THEN
           !
           ! ... final scf calculation with G-vectors for final cell
           !
           CALL reset_gvectors ( )
           !
        ELSE IF ( ions_status == 2 ) THEN
           !
           ! ... check whether nonzero magnetization is real
           !
           CALL reset_magn ( )
           !
        ELSE
           !
           ! ... update the wavefunctions, charge density, potential
           ! ... update_pot initializes structure factor array as well
           !
           CALL update_pot()
           !
           ! ... re-initialize atomic position-dependent quantities
           !
           CALL hinit1()
           !
        END IF
        !
        !
     END IF
     ! ... Reset convergence threshold of iterative diagonalization for
     ! ... the first scf iteration of each ionic step (after the first)
     !
     ethr = 1.0D-6
     !
  END DO main_loop
  !
  ! ... save final data file
  !
  CALL qexsd_set_status(exit_status)
  CALL punch('all')
  !
  CALL qmmm_shutdown()
  !
  IF ( .NOT. conv_ions )  exit_status =  3
  RETURN
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,&
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
  !
END SUBROUTINE run_pwscf

SUBROUTINE reset_gvectors ( )
  !
  ! ... Variable-cell optimization: once convergence is achieved, 
  ! ... make a final calculation with G-vectors and plane waves
  ! ... calculated for the final cell (may differ from the current
  ! ... result, using G_vectors and PWs for the starting cell)
  !
  USE io_global,  ONLY : stdout
  USE cellmd,     ONLY : lmovecell
  USE basis,      ONLY : starting_wfc, starting_pot
  USE fft_base,   ONLY : dfftp
  USE fft_base,   ONLY : dffts
  USE control_flags, ONLY : lbfgs, lmd
  IMPLICIT NONE
  !
  WRITE( UNIT = stdout, FMT = 9110 )
  WRITE( UNIT = stdout, FMT = 9120 )
  !
  ! ... prepare for a new scf, restarted from scratch, not from previous
  ! ... data (dimensions and file lengths will be different in general)
  !
  ! ... get magnetic moments from previous run before charge is deleted
  !
  CALL reset_starting_magnetization ( )
  !
  CALL clean_pw( .FALSE. )
  CALL close_files(.TRUE.)
  lmovecell=.FALSE.
  lbfgs=.FALSE.
  lmd=.FALSE.
  if (trim(starting_wfc) == 'file') starting_wfc = 'atomic+random'
  starting_pot='atomic'
  !
  ! ... re-set FFT grids
  !
  dfftp%nr1=0; dfftp%nr2=0; dfftp%nr3=0
  dffts%nr1=0; dffts%nr2=0; dffts%nr3=0
  !
  CALL init_run()
  !
9110 FORMAT( /5X,'A final scf calculation at the relaxed structure.' )
9120 FORMAT(  5X,'The G-vectors are recalculated for the final unit cell'/ &
              5X,'Results may differ from those at the preceding step.' )
  !
END SUBROUTINE reset_gvectors

SUBROUTINE reset_magn ( )
  !
  ! ... lsda optimization :  a final configuration with zero 
  ! ... absolute magnetization has been found and we check 
  ! ... if it is really the minimum energy structure by 
  ! ... performing a new scf iteration without any "electronic" history
  !
  USE io_global,  ONLY : stdout
  USE dfunct,     ONLY : newd
  IMPLICIT NONE
  !
  WRITE( UNIT = stdout, FMT = 9010 )
  WRITE( UNIT = stdout, FMT = 9020 )
  !
  ! ... re-initialize the potential (no need to re-initialize wavefunctions)
  !
  CALL potinit()
  CALL newd()

9010 FORMAT( /5X,'lsda relaxation :  a final configuration with zero', &
           & /5X,'                   absolute magnetization has been found' )
9020 FORMAT( /5X,'the program is checking if it is really ', &
           &     'the minimum energy structure',             &
           & /5X,'by performing a new scf iteration ',       & 
           &     'without any "electronic" history' )               
  !
END SUBROUTINE reset_magn
!
SUBROUTINE reset_starting_magnetization ( ) 
  !
  ! On input,  the scf charge density is needed
  ! On output, new values for starting_magnetization, angle1, angle2
  ! estimated from atomic magnetic moments - to be used in last step
  !
  USE kinds,     ONLY : dp
  USE constants, ONLY : pi
  USE ions_base, ONLY : nsp, ityp, nat
  USE lsda_mod,  ONLY : nspin, starting_magnetization
  USE scf,       ONLY : rho
  USE spin_orb,  ONLY : domag
  USE noncollin_module, ONLY : noncolin, angle1, angle2
  !
  IMPLICIT NONE
  INTEGER :: i, nt, iat
  REAL(dp):: norm_tot, norm_xy, theta, phi
  REAL (DP), ALLOCATABLE :: r_loc(:), m_loc(:,:)
  !
  IF ( (noncolin .AND. domag) .OR. nspin==2) THEN
     ALLOCATE ( r_loc(nat), m_loc(nspin-1,nat) )
     CALL get_locals(r_loc,m_loc,rho%of_r)
  ELSE
     RETURN
  END IF
  DO i = 1, nsp
     !
     starting_magnetization(i) = 0.0_DP
     angle1(i) = 0.0_DP
     angle2(i) = 0.0_DP
     nt = 0
     DO iat = 1, nat
        IF (ityp(iat) == i) THEN
           nt = nt + 1
           IF (noncolin) THEN
              norm_tot= sqrt(m_loc(1,iat)**2+m_loc(2,iat)**2+m_loc(3,iat)**2)
              norm_xy = sqrt(m_loc(1,iat)**2+m_loc(2,iat)**2)
              IF (norm_tot > 1.d-10) THEN
                 theta = acos(m_loc(3,iat)/norm_tot)
                 IF (norm_xy > 1.d-10) THEN
                    phi = acos(m_loc(1,iat)/norm_xy)
                    IF (m_loc(2,iat).lt.0.d0) phi = - phi
                 ELSE
                    phi = 2.d0*pi
                 END IF
              ELSE
                 theta = 2.d0*pi
                 phi = 2.d0*pi
              END IF
              angle1(i) = angle1(i) + theta
              angle2(i) = angle2(i) + phi
              starting_magnetization(i) = starting_magnetization(i) + &
                   norm_tot/r_loc(iat)
           ELSE
              starting_magnetization(i) = starting_magnetization(i) + &
                   m_loc(1,iat)/r_loc(iat)
           END IF
        END IF
     END DO
     IF ( nt > 0 ) THEN
        starting_magnetization(i) = starting_magnetization(i) / DBLE(nt)
        angle1(i) = angle1(i) / DBLE(nt)
        angle2(i) = angle2(i) / DBLE(nt)
     END IF
  END DO
  DEALLOCATE ( r_loc, m_loc )

END SUBROUTINE reset_starting_magnetization
