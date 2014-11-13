!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE move_ions()
  !----------------------------------------------------------------------------
  !
  ! ... This routine moves the ions according to the requested scheme:
  !
  ! ... lbfgs               bfgs minimizations
  ! ... lmd                 molecular dynamics ( verlet of vcsmd )
  ! ... lmd + lconstrain    constrained molecular dynamics,
  !
  ! ... coefficients for potential and wavefunctions extrapolation are
  ! ... also computed here
  !
  USE constants,              ONLY : e2, eps8, ry_kbar
  USE io_global,              ONLY : stdout
  USE io_files,               ONLY : tmp_dir, iunupdate, seqopn
  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : alat, at, bg, omega, cell_force, fix_volume, fix_area
  USE cellmd,                 ONLY : omega_old, at_old, press, lmovecell, calc
  USE ions_base,              ONLY : nat, ityp, tau, if_pos
  USE fft_base,               ONLY : dfftp
  USE fft_base,               ONLY : dffts
  USE grid_subroutines,       ONLY : realspace_grids_init
  USE gvect,                  ONLY : gcutm
  USE gvecs,                  ONLY : gcutms
  USE grid_subroutines,       ONLY : realspace_grids_init
  USE symm_base,              ONLY : checkallsym
  USE ener,                   ONLY : etot
  USE force_mod,              ONLY : force, sigma
  USE control_flags,          ONLY : istep, nstep, upscale, lbfgs, ldamped, &
                                     lconstrain, conv_ions, use_SMC, &
                                     lmd, llang, history, tr2
  USE basis,                  ONLY : starting_wfc
  USE relax,                  ONLY : epse, epsf, epsp, starting_scf_threshold
  USE lsda_mod,               ONLY : lsda, absmag
  USE mp_images,              ONLY : intra_image_comm
  USE io_global,              ONLY : ionode_id, ionode
  USE mp,                     ONLY : mp_bcast
  USE bfgs_module,            ONLY : bfgs, terminate_bfgs
  USE basic_algebra_routines, ONLY : norm
  USE dynamics_module,        ONLY : verlet, langevin_md, proj_verlet
  USE dynamics_module,        ONLY : smart_MC
  USE dfunct,                 only : newd
  !
  IMPLICIT NONE
  !
  LOGICAL, SAVE         :: lcheck_mag = .TRUE., &
                           restart_with_starting_magnetiz = .FALSE., &
                           lcheck_cell= .TRUE., &
                           final_cell_calculation=.FALSE.
  REAL(DP), ALLOCATABLE :: tauold(:,:,:)
  REAL(DP)              :: energy_error, gradient_error, cell_error
  LOGICAL               :: step_accepted, exst
  REAL(DP), ALLOCATABLE :: pos(:), grad(:)
  REAL(DP)              :: h(3,3), fcell(3,3)=0.d0, epsp1
  INTEGER,  ALLOCATABLE :: fixion(:)
  real(dp) :: tr
  !
  IF (use_SMC) CALL smart_MC()  ! for smart monte carlo method
  !
  ! ... only one node does the calculation in the parallel case
  !
  IF ( ionode ) THEN
     !
     conv_ions = .FALSE.
     !
     ALLOCATE( tauold( 3, nat, 3 ) )
     !
     ! ... the file containing old positions is opened 
     ! ... ( needed for extrapolation )
     !
     CALL seqopn( iunupdate, 'update', 'FORMATTED', exst ) 
     !
     IF ( exst ) THEN
        !
        READ( UNIT = iunupdate, FMT = * ) history
        READ( UNIT = iunupdate, FMT = * ) tauold
        !
     ELSE
        !
        history = 0
        tauold  = 0.D0
        !
        WRITE( UNIT = iunupdate, FMT = * ) history
        WRITE( UNIT = iunupdate, FMT = * ) tauold
        !
     END IF
     !
     CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
     !
     ! ... save the previous two steps ( a total of three steps is saved )
     !
     tauold(:,:,3) = tauold(:,:,2)
     tauold(:,:,2) = tauold(:,:,1)
     tauold(:,:,1) = tau(:,:)
     !
     ! ... history is updated (a new ionic step has been done)
     !
     history = MIN( 3, ( history + 1 ) )
     !
     ! ... old positions are written on file
     !
     CALL seqopn( iunupdate, 'update', 'FORMATTED', exst ) 
     !
     WRITE( UNIT = iunupdate, FMT = * ) history
     WRITE( UNIT = iunupdate, FMT = * ) tauold
     !
     CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
     !  
     DEALLOCATE( tauold )
     !
     ! ... do the minimization / dynamics step
     !
     IF ( lmovecell .AND. lconstrain ) THEN
        !
        IF ( lbfgs) CALL errore('move_ions', &
            & 'variable-cell bfgs and constraints not implemented yet', 1 )
        WRITE(stdout, '(5x,"-------------------------------------------")')
        WRITE(stdout, '(5x,"NEW FEATURE: constraints with variable cell")')
        WRITE(stdout, '(5x,"-------------------------------------------")')
        !
     END IF
     !
     ! ... BFGS algorithm is used to minimize ionic configuration
     !
     bfgs_minimization : &
     IF ( lbfgs ) THEN
        !
        ! ... the bfgs procedure is used
        !  
        ALLOCATE( pos( 3*nat ), grad( 3*nat ), fixion( 3*nat ) )
        !
        h = at * alat
        !
        pos    =   RESHAPE( tau,    (/ 3 * nat /) )
        CALL cryst_to_cart( nat, pos, bg, -1 )
        grad   = - RESHAPE( force,  (/ 3 * nat /) ) * alat
        CALL cryst_to_cart( nat, grad, at, -1 )
        fixion =   RESHAPE( if_pos, (/ 3 * nat /) )
        !
        IF ( lmovecell ) THEN
           at_old = at
           omega_old = omega
           etot = etot + press * omega
           CALL cell_force( fcell, - transpose(bg)/alat, sigma, omega, press )
           epsp1 = epsp / ry_kbar
        END IF
        !
        CALL bfgs( pos, h, etot, grad, fcell, fixion, tmp_dir, stdout, epse,&
                   epsf, epsp1,  energy_error, gradient_error, cell_error,  &
                   istep, nstep, step_accepted, conv_ions, lmovecell )
        !
        IF ( lmovecell ) THEN
           ! changes needed only if cell moves
           if (fix_volume) call impose_deviatoric_strain(alat*at, h)
           if (fix_area)   call impose_deviatoric_strain_2d(alat*at, h)
           at = h /alat  
           CALL recips( at(1,1),at(1,2),at(1,3), bg(1,1),bg(1,2),bg(1,3) )
           CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )
        END IF
        !
        CALL cryst_to_cart( nat, pos, at, 1 )
        tau    =   RESHAPE( pos, (/ 3 , nat /) )
        CALL cryst_to_cart( nat, grad, bg, 1 )
        force = - RESHAPE( grad, (/ 3, nat /) )
        !
        IF ( conv_ions ) THEN
           !
           IF ( ( lsda .AND. ( absmag < eps8 ) .AND. lcheck_mag ) ) THEN
              !
              ! ... lsda relaxation :  a final configuration with zero 
              ! ...                    absolute magnetization has been found.
              !                        A check on this configuration is needed
              restart_with_starting_magnetiz = .true.
              ! 
           ELSE IF (lmovecell.and.lcheck_cell) THEN
              !
              !  After the cell relaxation we make a final calculation
              !  with the correct g vectors corresponding to the relaxed
              !  cell.
              !
              final_cell_calculation=.TRUE.
              CALL terminate_bfgs ( etot, epse, epsf, epsp, lmovecell, &
                                    stdout, tmp_dir )
              !
           ELSE
              !
              CALL terminate_bfgs ( etot, epse, epsf, epsp, lmovecell, &
                                    stdout, tmp_dir )
              !
           END IF
           !
        ELSE
           !
           ! ... if a new bfgs step is done, new threshold is computed
           !
           IF ( step_accepted ) THEN
              !
              tr2  = starting_scf_threshold * &
                     MIN( 1.D0, ( energy_error / ( epse * upscale ) ), &
                                ( gradient_error / ( epsf * upscale ) ) )
              tr2  = MAX( ( starting_scf_threshold / upscale ), tr2 ) 
              !
           END IF
           !
           IF ( tr2 > 1.D-10 ) THEN
              WRITE( stdout, &
                     '(5X,"new conv_thr",T30,"= ",0PF18.10," Ry",/)' ) tr2
           ELSE
              WRITE( stdout, &
                     '(5X,"new conv_thr",T30,"= ",1PE18.1 ," Ry",/)' ) tr2
           END IF
           !
           ! ... the logical flag lcheck_mag is set again to .TRUE. (needed if 
           ! ... a new configuration with zero absolute magnetization is 
           ! ... identified in the following steps of the relaxation)
           !
           lcheck_mag = .TRUE.
           IF (lmovecell) lcheck_cell = .TRUE.
           !
        END IF
        !
        CALL output_tau( lmovecell, conv_ions )
        !
        DEALLOCATE( pos, grad, fixion )
        !
     END IF bfgs_minimization
     !
     ! ... molecular dynamics schemes are used
     !
     IF ( lmd ) THEN
        !
        IF ( calc == ' ' ) THEN
           !
           ! ... dynamics algorithms
           !
           IF ( ldamped ) THEN
              !
              CALL proj_verlet()
              !
           ELSE IF ( llang ) THEN
              !
              CALL langevin_md()
              !
           ELSE
              !
              CALL verlet()
              !
           END IF
           !
        ELSE IF ( calc /= ' ' ) THEN
           !
           ! ... variable cell shape md
           !
           CALL vcsmd()
           !
        END IF
        !
     END IF
     !
     ! ... before leaving check that the new positions still transform
     ! ... according to the symmetry of the system.
     !
     CALL checkallsym( nat, tau, ityp, dfftp%nr1, dfftp%nr2, dfftp%nr3 )
     !
  END IF

  CALL mp_bcast(restart_with_starting_magnetiz,ionode_id,intra_image_comm)
  CALL mp_bcast(final_cell_calculation,ionode_id,intra_image_comm)
  !
  IF ( final_cell_calculation ) THEN
     ! 
     ! ... Variable-cell optimization: once convergence is achieved, 
     ! ... make a final calculation with G-vectors and plane waves
     ! ... calculated for the final cell (may differ from the current
     ! ... result, using G_vectors and PWs for the starting cell)
     !
     WRITE( UNIT = stdout, FMT = 9110 )
     WRITE( UNIT = stdout, FMT = 9120 )
     !
     ! ... prepare for a new run, restarted from scratch, not from previous
     ! ... data (dimensions and file lengths will be different in general)
     !
     CALL clean_pw( .FALSE. )
     CALL close_files(.TRUE.)
     lmovecell=.FALSE.
     lcheck_cell=.FALSE.
     final_cell_calculation=.FALSE.
     lbfgs=.FALSE.
     lmd=.FALSE.
     lcheck_mag = .FALSE.
     restart_with_starting_magnetiz = .FALSE.
     if (trim(starting_wfc) == 'file') starting_wfc = 'atomic+random'
     ! ... conv_ions is set to .FALSE. to perform a final scf cycle
     conv_ions = .FALSE.
     ! ... allow re-calculation of FFT grid
     !
     dfftp%nr1=0; dfftp%nr2=0; dfftp%nr3=0; dffts%nr1=0; dffts%nr2=0; dffts%nr3=0
     CALL realspace_grids_init (dfftp, dffts,at, bg, gcutm, gcutms )
     CALL init_run()
     !
  ELSE IF (restart_with_starting_magnetiz) THEN
     !
     ! ... lsda optimization :  a final configuration with zero 
     ! ... absolute magnetization has been found and we check 
     ! ... if it is really the minimum energy structure by 
     ! ... performing a new scf iteration without any "electronic" history
     !
     WRITE( UNIT = stdout, FMT = 9010 )
     WRITE( UNIT = stdout, FMT = 9020 )
     !
     lcheck_mag = .FALSE.
     restart_with_starting_magnetiz = .FALSE.
     ! ... conv_ions is set to .FALSE. to perform a final scf cycle
     conv_ions = .FALSE.
     !
     ! ... re-initialize the potential and wavefunctions
     !
     CALL potinit()
     CALL newd()
     CALL wfcinit()
     !
  END IF
  !
  ! ... broadcast calculated quantities to all nodes
  !
  CALL mp_bcast( istep,     ionode_id, intra_image_comm )
  CALL mp_bcast( tau,       ionode_id, intra_image_comm )
  CALL mp_bcast( force,     ionode_id, intra_image_comm )
  CALL mp_bcast( tr2,       ionode_id, intra_image_comm )
  CALL mp_bcast( conv_ions, ionode_id, intra_image_comm )
  CALL mp_bcast( history,   ionode_id, intra_image_comm )
  !
  IF ( lmovecell ) THEN
     !
     CALL mp_bcast( at,        ionode_id, intra_image_comm )
     CALL mp_bcast( at_old,    ionode_id, intra_image_comm )
     CALL mp_bcast( omega,     ionode_id, intra_image_comm )
     CALL mp_bcast( omega_old, ionode_id, intra_image_comm )
     CALL mp_bcast( bg,        ionode_id, intra_image_comm )
     !
  END IF
  !
  RETURN

9010 FORMAT( /5X,'lsda relaxation :  a final configuration with zero', &
           & /5X,'                   absolute magnetization has been found' )
9020 FORMAT( /5X,'the program is checking if it is really ', &
           &     'the minimum energy structure',             &
           & /5X,'by performing a new scf iteration ',       & 
           &     'without any "electronic" history' )               
  !
9110 FORMAT( /5X,'A final scf calculation at the relaxed structure.' )
9120 FORMAT(  5X,'The G-vectors are recalculated for the final unit cell'/ &
              5X,'Results may differ from those at the preceding step.' )
  !
END SUBROUTINE move_ions
  !
