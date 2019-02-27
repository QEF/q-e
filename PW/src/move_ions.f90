!
! Copyright (C) 2002-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE move_ions ( idone, ions_status )
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
                                     fix_volume, fix_area, ibrav, enforce_ibrav
  USE cellmd,                 ONLY : omega_old, at_old, press, lmovecell, calc
  USE ions_base,              ONLY : nat, ityp, zv, tau, if_pos
  USE symm_base,              ONLY : checkallsym
  USE ener,                   ONLY : etot, ef
  USE force_mod,              ONLY : force, sigma
  USE control_flags,          ONLY : istep, nstep, upscale, lbfgs, &
                                     lconstrain, lmd, tr2
  USE relax,                  ONLY : epse, epsf, epsp, starting_scf_threshold
  USE lsda_mod,               ONLY : lsda, absmag
  USE mp_images,              ONLY : intra_image_comm
  USE io_global,              ONLY : ionode_id, ionode
  USE mp,                     ONLY : mp_bcast
  USE bfgs_module,            ONLY : bfgs, terminate_bfgs
  USE basic_algebra_routines, ONLY : norm
  USE dynamics_module,        ONLY : verlet, terminate_verlet, proj_verlet
  USE dynamics_module,        ONLY : smart_MC, langevin_md
  USE fcp                ,    ONLY : fcp_verlet, fcp_line_minimisation, &
                                     fcp_mdiis_update, fcp_mdiis_end
  USE fcp_variables,          ONLY : lfcpopt, lfcpdyn, fcp_mu, &
                                     fcp_relax, fcp_relax_crit
  USE klist,                  ONLY : nelec
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)   :: idone
  !! idone: see run_pwscf
  INTEGER,  INTENT(INOUT):: ions_status
  !! ions_status: see run_pwscf
  !
  REAL(DP)              :: energy_error, gradient_error, cell_error
  LOGICAL               :: step_accepted, exst
  REAL(DP), ALLOCATABLE :: pos(:), grad(:)
  REAL(DP)              :: h(3,3), fcell(3,3)=0.d0, epsp1
  INTEGER,  ALLOCATABLE :: fixion(:)
  LOGICAL               :: conv_fcp, conv_ions
  !
  ! ... only one node does the calculation in the parallel case
  !
  IF ( ionode ) THEN
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
     bfgs_minimization : &
     IF ( lbfgs ) THEN
        !
        ! ... BFGS algorithm is used to minimize ionic configuration
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
           IF ( TRIM(fcp_relax) == 'lm' ) THEN
              CALL fcp_line_minimisation( conv_fcp )
           ELSE IF ( TRIM(fcp_relax) == 'mdiis' ) THEN
              CALL fcp_mdiis_update( conv_fcp )
           END IF
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
           IF(enforce_ibrav) CALL remake_cell(ibrav, alat, at(1,1),at(1,2),at(1,3))  
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
           IF ( ions_status == 3 ) THEN
              !
              IF ( lsda .AND. absmag < eps6 ) THEN
                 !
                 ! ... a final configuration with zero absolute magnetization
                 ! ... has been found - do check with nonzero magnetization
                 !
                 ions_status = 2
                 !
              ELSE IF ( lmovecell ) THEN
                 !
                 ! ... Variable-cell relaxation converged with starting cell
                 ! ... Do final calculation with G-vectors for relaxed cell
                 !
                 ions_status = 1
                 !
              ELSE
                 !
                 ! ... Fixed-cell relaxation converged, prepare to exit
                 !
                 ions_status = 0
                 !
              END IF
              !
           ELSE IF ( ions_status == 2 ) THEN
              !
              ! ... check with nonzero magnetization succeeded, see above
              !
              IF ( lmovecell ) THEN
                 ions_status = 1
              ELSE
                 ions_status = 0
              END IF
              !
           ELSE IF ( ions_status == 1 ) THEN
              !
              ions_status = 0
              !
           END IF
           !
           IF ( ions_status < 2 ) THEN
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
              IF ( TRIM(fcp_relax) == 'mdiis' ) THEN
                 CALL fcp_mdiis_end()
              END IF
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
        END IF
        !
        CALL output_tau( lmovecell, conv_ions )
        !
        DEALLOCATE( pos, grad, fixion )
        !
     END IF bfgs_minimization
     !
     IF ( lmd ) THEN
        !
        conv_ions = .FALSE.
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
              IF ( TRIM(fcp_relax) == 'lm' ) THEN
                 CALL fcp_line_minimisation( conv_fcp )
              ELSE IF ( TRIM(fcp_relax) == 'mdiis' ) THEN
                 CALL fcp_mdiis_update( conv_fcp )
              END IF
              IF ( .not. conv_fcp .and. idone < nstep ) conv_ions = .FALSE.
              !
              ! ... FCP output
              !
              IF ( conv_ions ) THEN
                 WRITE( stdout, '(5X,"FCP : converged ", &
                      & "( criteria force < ", ES8.1," )")') fcp_relax_crit
                 WRITE( stdout, '(5X,"FCP : final tot_charge =",F12.6,/)') &
                      SUM( zv(ityp(1:nat)) ) - nelec
                 IF ( TRIM(fcp_relax) == 'mdiis' ) THEN
                    CALL fcp_mdiis_end()
                 END IF
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
        IF ( conv_ions ) ions_status  = 0
        !
     END IF
     !
     ! ... before leaving check that the new positions still transform
     ! ... according to the symmetry of the system.
     !
     CALL checkallsym( nat, tau, ityp)
     !
  END IF
  !
  
  CALL mp_bcast( ions_status, ionode_id, intra_image_comm )
  IF ( lfcpopt .or. lfcpdyn ) CALL mp_bcast(nelec,ionode_id,intra_image_comm)
  !
  !
  ! ... broadcast calculated quantities to all nodes
  !
  CALL mp_bcast( istep,     ionode_id, intra_image_comm )
  CALL mp_bcast( tau,       ionode_id, intra_image_comm )
  CALL mp_bcast( force,     ionode_id, intra_image_comm )
  CALL mp_bcast( tr2,       ionode_id, intra_image_comm )
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
  !
END SUBROUTINE move_ions
