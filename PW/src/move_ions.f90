!
! Copyright (C) 2002-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE move_ions ( idone )
  !----------------------------------------------------------------------------
  !
  ! ... Perform an ionic step, according to the requested scheme:
  ! ...    lbfgs               bfgs minimizations
  ! ...    lmd                 molecular dynamics ( all kinds )
  ! ... Additional variables affecting the calculation:
  ! ...    lmovecell           Variable-cell calculation
  ! ...    calc                Type of MD
  ! ...    lconstrain          constrained MD
  ! ..  "idone" is the counter on ionic moves, "nstep" their total number 
  ! ... "istep" contains the number of all steps including previous runs
  ! ... Coefficients for potential and wavefunctions extrapolation are
  ! ... no longer computed here but in update_pot
  !
  USE constants,              ONLY : e2, eps6, ry_kbar
  USE io_global,              ONLY : stdout
  USE io_files,               ONLY : tmp_dir
  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : alat, at, bg, omega, cell_force, &
                                     fix_volume, fix_area
  USE cellmd,                 ONLY : omega_old, at_old, press, lmovecell, calc
  USE ions_base,              ONLY : nat, ityp, zv, tau, if_pos
  USE fft_base,               ONLY : dfftp
  USE fft_base,               ONLY : dffts
  USE fft_types,              ONLY : fft_type_allocate
  USE gvect,                  ONLY : gcutm
  USE gvecs,                  ONLY : gcutms
  USE symm_base,              ONLY : checkallsym
  USE ener,                   ONLY : etot, ef
  USE force_mod,              ONLY : force, sigma
  USE control_flags,          ONLY : istep, nstep, upscale, lbfgs, &
                                     lconstrain, conv_ions, lmd, tr2
  USE basis,                  ONLY : starting_wfc
  USE relax,                  ONLY : epse, epsf, epsp, starting_scf_threshold
  USE lsda_mod,               ONLY : lsda, absmag
  USE mp_images,              ONLY : intra_image_comm
  USE mp_bands,               ONLY : intra_bgrp_comm
  USE io_global,              ONLY : ionode_id, ionode
  USE mp,                     ONLY : mp_bcast
  USE bfgs_module,            ONLY : bfgs, terminate_bfgs
  USE basic_algebra_routines, ONLY : norm
  USE dynamics_module,        ONLY : verlet, terminate_verlet, proj_verlet
  USE dynamics_module,        ONLY : smart_MC, langevin_md
  USE fcp                ,    ONLY : fcp_verlet, fcp_line_minimisation
  USE fcp_variables,          ONLY : lfcpopt, lfcpdyn, fcp_mu, &
                                     fcp_relax_crit
  USE klist,                  ONLY : nelec
  USE dfunct,                 only : newd
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: idone
  !
  LOGICAL, SAVE         :: lcheck_mag = .TRUE., &
                           restart_with_starting_magnetiz = .FALSE., &
                           lcheck_cell= .TRUE., &
                           final_cell_calculation=.FALSE.
  REAL(DP)              :: energy_error, gradient_error, cell_error
  LOGICAL               :: step_accepted, exst
  REAL(DP), ALLOCATABLE :: pos(:), grad(:)
  REAL(DP)              :: h(3,3), fcell(3,3)=0.d0, epsp1
  INTEGER,  ALLOCATABLE :: fixion(:)
  LOGICAL               :: conv_fcp
  !
  ! ... only one node does the calculation in the parallel case
  !
  IF ( ionode ) THEN
     !
     conv_ions = .FALSE.
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
                   lmovecell, step_accepted, conv_ions, istep )
        !
        ! ... relax for FCP
        !
        IF ( lfcpopt ) THEN
           CALL fcp_line_minimisation( conv_fcp )
           IF ( .not. conv_fcp .and. idone < nstep ) THEN
             conv_ions = .FALSE.
           END IF
        END IF
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
           IF ( ( lsda .AND. ( absmag < eps6 ) .AND. lcheck_mag ) ) THEN
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
           ! ... FCP output
           !
           IF ( lfcpopt ) THEN
             WRITE( stdout, '(/,5X, "FCP Optimisation : converged ", &
               & "( criteria force < ",ES8.1," )")') fcp_relax_crit
             WRITE( stdout, '(5X,"FCP Optimisation : tot_charge =",F12.6,/)') &
               SUM( zv(ityp(1:nat)) ) - nelec
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
        ! ... fixed-cell molecular dynamics algorithms first:
        ! ... projected Verlet, Langevin, Verlet
        !
        IF ( calc == 'vm' ) THEN
           !
           CALL proj_verlet( conv_ions )
           !
           ! ... relax for FCP
           !
           IF ( lfcpopt ) THEN
              CALL fcp_line_minimisation( conv_fcp )
              IF ( .not. conv_fcp .and. idone < nstep ) conv_ions = .FALSE.
              !
              ! ... FCP output
              !
              IF ( conv_ions ) THEN
                 WRITE( stdout, '(5X,"FCP : converged ", &
                      & "( criteria force < ", ES8.1," )")') fcp_relax_crit
                 WRITE( stdout, '(5X,"FCP : final tot_charge =",F12.6,/)') &
                      SUM( zv(ityp(1:nat)) ) - nelec
              END IF
           END IF
           IF ( .NOT. conv_ions .AND. idone >= nstep ) THEN
              WRITE( UNIT = stdout, FMT =  &
                   '(/,5X,"The maximum number of steps has been reached.")' )
              WRITE( UNIT = stdout, &
                   FMT = '(/,5X,"End of molecular dynamics calculation")' )
           END IF
           !
        ELSE IF ( calc(1:1) == 'l' ) THEN
           !
           ! ... for smart monte carlo method
           !
           IF ( calc(2:2) == 's' ) CALL smart_MC()
           !
           CALL langevin_md()
           !
           ! ... dynamics for FCP
           !
           IF ( lfcpdyn ) CALL fcp_verlet()
           IF ( idone >= nstep ) THEN
              WRITE( UNIT = stdout, FMT =  &
                   '(/,5X,"The maximum number of steps has been reached.")' )
              WRITE( UNIT = stdout, &
                   FMT = '(/,5X,"End of molecular dynamics calculation")' )
              conv_ions = .true.
           END IF
           !
        ELSE IF ( calc == 'vd' ) THEN
           !
           CALL verlet()
           !
           ! ... dynamics for FCP
           !
           IF ( lfcpdyn ) CALL fcp_verlet()
           IF ( idone >= nstep) THEN
              CALL terminate_verlet()
              conv_ions = .true.
           END IF
           !
        ELSE
           !
           ! ... variable cell shape md
           !
           CALL vcsmd( conv_ions )
           !
           ! ... after nstep, set conv_ions to T for MD, to F for damped MD
           !
           IF ( .NOT.conv_ions .AND. idone >= nstep ) THEN
              WRITE( UNIT = stdout, FMT = '(/,5X,"Maximum number of ", &
             &    "iterations reached, stopping")' )
              conv_ions = ( calc(2:2) == 'd' )
           END IF
           !
        END IF
        !
     END IF
     !
     ! ... before leaving check that the new positions still transform
     ! ... according to the symmetry of the system.
     !
     CALL checkallsym( nat, tau, ityp)
     !
  END IF

  CALL mp_bcast(restart_with_starting_magnetiz,ionode_id,intra_image_comm)
  CALL mp_bcast(final_cell_calculation,ionode_id,intra_image_comm)
  IF ( lfcpopt .or. lfcpdyn ) CALL mp_bcast(nelec,ionode_id,intra_image_comm)
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
     ! ... get magnetic moments from previous run before charge is deleted
     !
     CALL reset_starting_magnetization ( )
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
     !
     ! ... re-set and re-calculate FFT grid 
     !
     dfftp%nr1=0; dfftp%nr2=0; dfftp%nr3=0
     CALL fft_type_allocate (dfftp, at, bg, gcutm, intra_bgrp_comm )
     dffts%nr1=0; dffts%nr2=0; dffts%nr3=0
     CALL fft_type_allocate (dffts, at, bg, gcutms, intra_bgrp_comm)
     !
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
     ! ... re-initialize the potential (no need to re-initialize wavefunctions)
     !
     CALL potinit()
     CALL newd()
     !!! CALL wfcinit()
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
     starting_magnetization(i) = starting_magnetization(i) / REAL(nt)
     angle1(i) = angle1(i) / REAL(nt)
     angle2(i) = angle2(i) / REAL(nt)
  END DO
  DEALLOCATE ( r_loc, m_loc )

END SUBROUTINE reset_starting_magnetization
