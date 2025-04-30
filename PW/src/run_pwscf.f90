!
! Copyright (C) 2013-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE run_pwscf( exit_status ) 
  !----------------------------------------------------------------------------
  !! Author: Paolo Giannozzi  
  !! License: GNU  
  !! Summary: Run an instance of the Plane Wave Self-Consistent Field code
  !
  !! Run an instance of the Plane Wave Self-Consistent Field code 
  !! MPI initialization and input data reading is performed in the 
  !! calling code - returns in exit_status the exit code for pw.x:
  !! * 0: completed successfully
  !! * 1: an error has occurred (value returned by the errore() routine)
  !! * 2-127: convergence error
  !!    * 2: scf convergence error
  !!    * 3: ion convergence error
  !! * 128-255: code exited due to specific trigger
  !!    * 255: exit due to user request, or signal trapped,
  !!          or time > max_seconds
  !!     (note: in the future, check_stop_now could also return a value
  !!     to specify the reason of exiting, and the value could be used
  !!     to return a different value for different reasons)
  !! @Note
  !! 20/04/23 Unless preprocessing flag __RETURN_EXIT_STATUS is set (see
  !! routine do_stop) pw.x no longer returns an exit status != 0  to the shell
  !! @endnote
  !
  !! @Note
  !! 10/01/17 Samuel Ponce: Add Ford documentation
  !! @endnote
  !!
  !
  USE kinds,                ONLY : DP
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE parameters,           ONLY : ntypx, npk
  USE upf_params,           ONLY : lmaxx
  USE cell_base,            ONLY : fix_volume, fix_area
  USE control_flags,        ONLY : conv_elec, gamma_only, ethr, lscf, treinit_gvecs
  USE control_flags,        ONLY : conv_ions, istep, nstep, restart, lmd, lbfgs,&
                                   lensemb, lforce=>tprnfor, tstress
  USE cellmd,               ONLY : lmovecell
  USE command_line_options, ONLY : command_line
  USE force_mod,            ONLY : sigma, force
  USE check_stop,           ONLY : check_stop_init, check_stop_now
  USE mp_images,            ONLY : intra_image_comm
  USE extrapolation,        ONLY : update_file, update_pot
  USE scf,                  ONLY : rho
  USE lsda_mod,             ONLY : nspin
  USE fft_base,             ONLY : dfftp
  USE qmmm,                 ONLY : qmmm_initialization, qmmm_shutdown, &
                                   qmmm_update_positions, qmmm_update_forces
  USE qexsd_module,         ONLY : qexsd_set_status
  USE xc_lib,               ONLY : xclib_dft_is, stop_exx, exx_is_active
  USE beef,                 ONLY : beef_energies
  USE ldaU,                 ONLY : lda_plus_u
  USE add_dmft_occ,         ONLY : dmft
  USE extffield,            ONLY : init_extffield, close_extffield
  USE input_parameters,     ONLY : nextffield
  !
  USE device_fbuff_m,             ONLY : dev_buf
  !
#if defined (__ENVIRON)
  USE plugin_flags,      ONLY : use_environ
  USE environ_pw_module, ONLY : is_ms_gcs, init_ms_gcs
#endif
#if defined (__OSCDFT)
  USE plugin_flags,      ONLY : use_oscdft
  USE oscdft_base,       ONLY : oscdft_ctx
  USE oscdft_functions,  ONLY : oscdft_run_pwscf
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(OUT) :: exit_status
  !! Gives the exit status at the end
  !
  LOGICAL, EXTERNAL :: matches
  ! checks if first string is contained in the second
  !
  ! ... local variables
  !
  INTEGER :: idone 
  ! counter of electronic + ionic steps done in this run
  INTEGER :: ions_status
  ! ions_status =  3  not yet converged
  ! ions_status =  2  converged, restart with nonzero magnetization
  ! ions_status =  1  converged, final step with current cell needed
  ! ions_status =  0  converged, exiting
  !
  LOGICAL :: optimizer_failed = .FALSE.
  !
  INTEGER :: ierr
  ! collect error codes
  !
  ions_status = 3
  exit_status = 0
  IF ( ionode ) WRITE( UNIT = stdout, FMT = 9010 ) ntypx, npk, lmaxx
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
     IF (ionode) CALL run_dist( exit_status )
     RETURN
  ENDIF
  !
  IF ( gamma_only ) WRITE( UNIT = stdout, &
     & FMT = '(/,5X,"gamma-point specific algorithms are used")' )
  !
  ! call to void routine for user defined / plugin patches initializations
  !
#if defined(__LEGACY_PLUGINS)
  CALL plugin_initialization()
#endif 
#if defined (__ENVIRON)
  IF (use_environ) THEN
     IF (is_ms_gcs()) CALL init_ms_gcs()
  END IF
#endif
  !
  CALL check_stop_init()
  !
  CALL setup()
  !
  CALL qmmm_update_positions()
  !
  ! ... dry run: code will stop here if called with exit file present
  ! ... useful for a quick and automated way to check input data
  !
  IF ( nstep == 0 .OR. check_stop_now() ) THEN
     CALL pre_init()
     CALL data_structure( gamma_only )
     CALL summary()
     CALL memory_report()
     exit_status = 255
     CALL qexsd_set_status( exit_status )
     CALL punch( 'config-init' )
     RETURN
  ENDIF
  !
  CALL init_run()
  !
  !  read external force fields parameters
  ! 
  IF ( nextffield > 0 .AND. ionode) THEN
     !
     CALL init_extffield( 'PW', nextffield )
     !
  END IF
  !
  IF ( check_stop_now() ) THEN
     exit_status = 255
     CALL qexsd_set_status( exit_status )
     CALL punch( 'config' )
     RETURN
  ENDIF
  !
  main_loop: DO idone = 1, nstep
     !
     ! ... electronic self-consistency or band structure calculation
     !
#if defined (__OSCDFT)
     IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) THEN
        CALL oscdft_run_pwscf(oscdft_ctx)
     ELSE
#endif
     IF ( .NOT. lscf) THEN
        CALL non_scf()
     ELSE
        CALL electrons()
     END IF
#if defined (__OSCDFT)
     END IF
#endif
     !
     ! ... code stopped by user or not converged
     !
     IF ( check_stop_now() .OR. .NOT. conv_elec ) THEN
        IF ( check_stop_now() ) THEN
            exit_status = 255
        ELSE
           IF (dmft) THEN
              exit_status =  131
           ELSE
              exit_status = 2
           ENDIF
        ENDIF
        CALL qexsd_set_status(exit_status)
        IF(exx_is_active()) then
          CALL punch( 'all' )
        ELSE
          CALL punch( 'config' )
        ENDIF
        RETURN
     ENDIF
     !
     ! ... file in CASINO format written here if required
     !
     IF ( lmd ) THEN
        CALL pw2casino( istep )
     ELSE
        CALL pw2casino( 0 )
     END IF
     !
     ! ... ionic section starts here
     !
     CALL start_clock( 'ions' ); !write(*,*)' start ions' ; FLUSH(6)
     conv_ions = .TRUE.
     !
     ! ... force calculation
     !
     IF ( lforce ) CALL forces()
     !
     ! ... stress calculation
     !
     IF ( tstress ) CALL stress( sigma )
     !
     IF ( lmd .OR. lbfgs ) THEN
        !
        ! ... add information on this ionic step to xml file
        !
        CALL add_qexsd_step( idone )
        !
        IF (fix_volume) CALL impose_deviatoric_stress( sigma )
        IF (fix_area)   CALL impose_deviatoric_stress_2d( sigma )
        !
        ! ... save data needed for potential and wavefunction extrapolation
        !
        CALL update_file()
        !
        ! ... ionic step (for molecular dynamics or optimization)
        !
        CALL move_ions ( idone, ions_status, optimizer_failed )
        conv_ions = ( ions_status == 0 ) .OR. &
                    ( ions_status == 1 .AND. treinit_gvecs )
        !
        IF ( xclib_dft_is('hybrid') )  CALL stop_exx()
        !
        ! ... save restart information for the new configuration
        !
        IF ( idone <= nstep .AND. .NOT. conv_ions ) THEN
            exit_status = 255
            CALL qexsd_set_status( exit_status )
            CALL punch( 'config-only' )
        END IF
        !
     END IF
     !
     CALL stop_clock( 'ions' ); !write(*,*)' stop ions' ; FLUSH(6)
     !
     ! ... send out forces to MM code in QM/MM run
     !
     CALL qmmm_update_forces( force, rho%of_r, nspin, dfftp )
     !
     ! ... exit condition (ionic convergence) is checked here
     !
     IF ( conv_ions .OR. optimizer_failed ) EXIT main_loop
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
           lbfgs=.FALSE.; lmd=.FALSE.
           WRITE( UNIT = stdout, FMT=9020 ) 
           !
           CALL reset_gvectors( )
           !
           ! ... read atomic occupations for DFT+U(+V)
           !
           IF ( lda_plus_u ) CALL read_ns()
           !
        ELSE IF ( ions_status == 2 ) THEN
           !
           ! ... check whether nonzero magnetization is real
           !
           CALL reset_magn()
           !
        ELSE
           !
           IF ( treinit_gvecs ) THEN
              !
              ! ... prepare for next step with freshly computed G vectors
              !
              IF ( lmovecell) CALL scale_h()
              CALL reset_gvectors ( )
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
        END IF
        !
     ENDIF
     ! ... Reset convergence threshold of iterative diagonalization for
     ! ... the first scf iteration of each ionic step (after the first)
     !
     ethr = 1.0D-6
     !
     CALL dev_buf%reinit( ierr )
     IF ( ierr .ne. 0 ) CALL infomsg( 'run_pwscf', 'Cannot reset GPU buffers! Some buffers still locked.' )
     !
  ENDDO main_loop
  !
  ! Set correct exit_status
  !
  IF ( .NOT. conv_ions .OR. optimizer_failed ) THEN
      exit_status =  3
  ELSE
      ! All good
      exit_status = 0
   END IF
  !
  ! ... save final data file
  !
  CALL qexsd_set_status( exit_status )
  IF ( lensemb ) CALL beef_energies( )
  CALL punch( 'all' )
  !
  CALL qmmm_shutdown()
  !
  RETURN
  !
9010 FORMAT( /,5X,'Current dimensions of program PWSCF are:', &
           & /,5X,'Max number of different atomic species (ntypx) = ',I2,&
           & /,5X,'Max number of k-points (npk) = ',I6,       &
           & /,5X,'Max angular momentum in pseudopotentials (lmaxx) = ',i2)
9020 FORMAT( /,5X,'Final scf calculation at the relaxed structure.', &
          &  /,5X,'The G-vectors are recalculated for the final unit cell', &
          &  /,5X,'Results may differ from those at the preceding step.' )
  !
END SUBROUTINE run_pwscf
!
!
!-------------------------------------------------------------
SUBROUTINE reset_gvectors( )
!-------------------------------------------------------------
  !
  !! Prepare a new scf calculation with newly recomputed grids,
  !! restarting from scratch, not from available data of previous
  !! steps (dimensions and file lengths will be different in general)
  !! Useful as a check of variable-cell optimization: 
  !! once convergence is achieved, compare the final energy with the
  !! energy computed with G-vectors and plane waves for the final cell
  !
  USE io_global,  ONLY : stdout
  USE starting_scf, ONLY : starting_wfc, starting_pot
  USE fft_base,   ONLY : dfftp
  USE fft_base,   ONLY : dffts
  USE xc_lib,     ONLY : xclib_dft_is
  ! 
  IMPLICIT NONE
  !
  ! ... get magnetic moments from previous run before charge is deleted
  !
  CALL reset_starting_magnetization()
  !
  ! ... clean everything (FIXME: clean only what has to be cleaned)
  !
  CALL clean_pw( .FALSE. )
  CALL close_files(.TRUE.)
  !
  IF (TRIM(starting_wfc) == 'file') starting_wfc = 'atomic+random'
  starting_pot='atomic'
  !
  ! ... re-set FFT grids and re-compute needed stuff (FIXME: which?)
  !
  dfftp%nr1=0; dfftp%nr2=0; dfftp%nr3=0
  dffts%nr1=0; dffts%nr2=0; dffts%nr3=0
  !
  CALL init_run()
  !
  ! ... re-set and re-initialize EXX-related stuff
  !
  IF ( xclib_dft_is('hybrid') ) CALL reset_exx( )
  !
END SUBROUTINE reset_gvectors
!
!
!-------------------------------------------------------------
SUBROUTINE reset_exx( )
!-------------------------------------------------------------
  USE fft_types,  ONLY : fft_type_deallocate 
  USE exx_base,   ONLY : exx_grid_init, exx_mp_init, exx_div_check, & 
                         coulomb_fac, coulomb_done 
  USE exx,        ONLY : dfftt, exx_fft_create, deallocate_exx 
  USE exx_band,   ONLY : igk_exx 
  ! 
  IMPLICIT NONE
  !
  ! ... re-set EXX-related stuff...
  !
  IF (ALLOCATED(coulomb_fac) ) DEALLOCATE( coulomb_fac, coulomb_done )
  CALL deallocate_exx( )
  IF (ALLOCATED(igk_exx)) DEALLOCATE(igk_exx) 
  dfftt%nr1=0; dfftt%nr2=0; dfftt%nr3=0 
  CALL fft_type_deallocate( dfftt ) ! FIXME: is this needed?
  !
  ! ... re-compute needed EXX-related stuff
  !
  CALL exx_grid_init( REINIT = .TRUE. )
  CALL exx_mp_init()
  CALL exx_fft_create()
  CALL exx_div_check()
  ! 
END SUBROUTINE reset_exx
!
!
!----------------------------------------------------------------
SUBROUTINE reset_magn()
  !----------------------------------------------------------------
  !! LSDA optimization: a final configuration with zero 
  !! absolute magnetization has been found and we check 
  !! if it is really the minimum energy structure by 
  !! performing a new scf iteration without any "electronic" history.
  !
  USE io_global,    ONLY : stdout
  USE dfunct,       ONLY : newd
  !
  IMPLICIT NONE
  !
  WRITE( UNIT = stdout, FMT = 9010 )
  WRITE( UNIT = stdout, FMT = 9020 )
  !
  ! ... re-initialize the potential (no need to re-initialize wavefunctions)
  !
  CALL potinit()
  CALL newd()
  !
9010 FORMAT( /5X,'lsda relaxation :  a final configuration with zero', &
           & /5X,'                   absolute magnetization has been found' )
9020 FORMAT( /5X,'the program is checking if it is really ', &
           &     'the minimum energy structure',             &
           & /5X,'by performing a new scf iteration ',       & 
           &     'without any "electronic" history' )               
  !
END SUBROUTINE reset_magn
!
!
!-------------------------------------------------------------------
SUBROUTINE reset_starting_magnetization() 
  !-------------------------------------------------------------------
  !! On input, the scf charge density is needed.  
  !! On output, new values for starting_magnetization, angle1, angle2
  !! estimated from atomic magnetic moments - to be used in last step.
  !
  USE kinds,              ONLY : DP
  USE constants,          ONLY : pi
  USE ions_base,          ONLY : nsp, ityp, nat
  USE lsda_mod,           ONLY : nspin, starting_magnetization
  USE scf,                ONLY : rho
  USE noncollin_module,   ONLY : noncolin, angle1, angle2, domag
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER  :: i, nt, iat
  ! loop counter on species
  ! number of atoms per species
  ! loop counter on atoms
  REAL(DP) :: norm_tot, norm_xy
  ! modulus of atomic magnetization
  ! xy-projection of atomic magnetization
  REAL(DP) :: theta, phi
  ! angle between magnetization and z-axis
  ! angle between xy-magnetization and x-axis
  REAL(DP), ALLOCATABLE :: r_loc(:)
  ! auxiliary array for density
  REAL(DP), ALLOCATABLE :: m_loc(:,:)
  ! auxiliary array for magnetization
  !
  IF ( (noncolin .AND. domag) .OR. nspin==2) THEN
     ALLOCATE( r_loc(nat), m_loc(nspin-1,nat) )
     CALL get_locals( r_loc,m_loc, rho%of_r )
  ELSE
     RETURN
  ENDIF
  !
  DO i = 1, nsp
     !
     starting_magnetization(i) = 0.0_DP
     angle1(i) = 0.0_DP
     angle2(i) = 0.0_DP
     nt = 0
     !
     DO iat = 1, nat
        IF (ityp(iat) == i) THEN
           nt = nt + 1
           IF (noncolin) THEN
              norm_tot = SQRT(m_loc(1,iat)**2+m_loc(2,iat)**2+m_loc(3,iat)**2)
              norm_xy  = SQRT(m_loc(1,iat)**2+m_loc(2,iat)**2)
              IF (norm_tot > 1.d-10) THEN
                 theta = ACOS(m_loc(3,iat)/norm_tot)
                 IF (norm_xy > 1.d-10) THEN
                    phi = ACOS(m_loc(1,iat)/norm_xy)
                    IF (m_loc(2,iat) < 0.d0) phi = - phi
                 ELSE
                    phi = 2.d0*pi
                 ENDIF
              ELSE
                 theta = 2.d0*pi
                 phi = 2.d0*pi
              ENDIF
              angle1(i) = angle1(i) + theta
              angle2(i) = angle2(i) + phi
              starting_magnetization(i) = starting_magnetization(i) + &
                                          norm_tot/r_loc(iat)
           ELSE
              starting_magnetization(i) = starting_magnetization(i) + &
                                          m_loc(1,iat)/r_loc(iat)
           ENDIF
        ENDIF
     ENDDO
     !
     IF ( nt > 0 ) THEN
        starting_magnetization(i) = starting_magnetization(i) / DBLE(nt)
        angle1(i) = angle1(i) / DBLE(nt)
        angle2(i) = angle2(i) / DBLE(nt)
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE( r_loc, m_loc )
  !
END SUBROUTINE reset_starting_magnetization
