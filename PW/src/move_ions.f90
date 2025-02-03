!
! Copyright (C) 2002-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE move_ions( idone, ions_status, optimizer_failed )
  !----------------------------------------------------------------------------
  !! Perform a ionic step, according to the requested scheme:
  !
  !! * lbfgs: bfgs minimizations
  !! * lmd: molecular dynamics ( all kinds )
  !
  !! Additional variables affecting the calculation:
  !
  !! * lmovecell: Variable-cell calculation
  !! * calc: type of MD
  !! * lconstrain: constrained MD
  !! * "idone" is the counter on ionic moves, "nstep" their total number 
  !! * "istep" contains the number of all steps including previous runs.
  !
  !! Coefficients for potential and wavefunctions extrapolation are
  !! no longer computed here but in update_pot.
  !
  USE constants,              ONLY : e2, eps4, eps6, ry_kbar
  USE control_flags,          ONLY : tnosep
  USE io_global,              ONLY : stdout
  USE io_files,               ONLY : tmp_dir, prefix
  USE kinds,                  ONLY : DP
  USE cell_base,              ONLY : alat, at, bg, omega, cell_force, &
                                     fix_volume, fix_area, ibrav, press, &
                                     iforceh, enforce_ibrav
  USE cellmd,                 ONLY : omega_old, at_old, lmovecell, calc
  USE ions_base,              ONLY : nat, ityp, zv, tau, if_pos
  USE symm_base,              ONLY : checkallsym
  USE ener,                   ONLY : etot, ef
  USE force_mod,              ONLY : force, sigma
  USE control_flags,          ONLY : istep, nstep, upscale, lbfgs, &
                                     lconstrain, lmd, tr2, iprint, tnosep
  USE relax,                  ONLY : epse, epsf, epsp, starting_scf_threshold
  USE lsda_mod,               ONLY : lsda, absmag
  USE mp_images,              ONLY : intra_image_comm
  USE io_global,              ONLY : ionode_id, ionode
  USE mp,                     ONLY : mp_bcast
  USE bfgs_module,            ONLY : bfgs, terminate_bfgs
  USE ions_nose,              ONLY : ions_nosevel, ions_noseupd, ions_nose_shiftvar, vnhp, xnhp0, xnhpp, xnhpm, nhpcl, nhpdim,&
                                     kbt, nhpbeg, nhpend, ekin2nhp, qnp, gkbt2nhp, ions_nose_energy, ions_nose_nrg
  USE basic_algebra_routines, ONLY : norm
  USE dynamics_module,        ONLY : verlet, terminate_verlet, proj_verlet, fire
  USE dynamics_module,        ONLY : smart_MC, langevin_md, dt, vel, elapsed_time
  USE dynamics_module,        ONLY : fire_nmin, fire_f_inc, fire_f_dec, &
                                     fire_alpha_init, fire_falpha, fire_dtmax, RyDt_to_HaDt
  USE dynamics_module,        ONLY : velocity_verlet
  USE klist,                  ONLY : nelec, tot_charge
  USE dfunct,                 only : newd
  USE fcp_module,             ONLY : lfcp, fcp_eps, fcp_mu, fcp_relax, &
                                     fcp_verlet, fcp_terminate, output_fcp
  USE rism_module,            ONLY : lrism, rism_new_conv_thr
  USE printout_base,          ONLY : printout_base_open, printout_base_close, &
                                     printout_cell, printout_pos, printout_stress
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)   :: idone
  !! idone: see run_pwscf
  INTEGER,  INTENT(INOUT):: ions_status
  !! ions_status: see run_pwscf
  LOGICAL,  INTENT(OUT):: optimizer_failed
  !
  ! ... local variables
  !
  REAL(DP)              :: energy_error, gradient_error, cell_error, fcp_error
  LOGICAL               :: step_accepted, exst
  REAL(DP), ALLOCATABLE :: pos(:), grad(:)
  REAL(DP)              :: h(3,3), fcell(3,3)=0.d0, epsp1, new_alat
  REAL(DP)              :: relec, felec, helec, capacitance, tot_charge_
  LOGICAL               :: conv_ions
  CHARACTER(LEN=320)    :: filebfgs
  INTEGER               :: iunit
  INTEGER               :: nose_cycle
  !
  optimizer_failed = .FALSE.
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
     ENDIF
     !
     bfgs_minimization : &
     IF ( lbfgs ) THEN
        !
        ! ... BFGS algorithm is used to minimize ionic configuration
        !
        ALLOCATE( pos( 3*nat ), grad( 3*nat ) )
        pos = 0.d0; grad=0.d0
        !
        h = at * alat
        !
        pos    =   RESHAPE( tau,    (/ 3 * nat /) )
        CALL cryst_to_cart( nat, pos, bg, -1 )
        grad   = - RESHAPE( force,  (/ 3 * nat /) ) * alat
        CALL cryst_to_cart( nat, grad, at, -1 )
        !
        IF ( lmovecell ) THEN
           at_old = at
           omega_old = omega
           etot = etot + press * omega
           CALL cell_force( fcell, - transpose(bg)/alat, sigma, omega, press )
           epsp1 = epsp / ry_kbar
        ENDIF
        !
        relec = 0.0_DP
        felec = 0.0_DP
        IF ( lfcp ) THEN
           relec = nelec
           felec = (ef - fcp_mu)
           CALL fcp_capacitance( capacitance, -1.0_DP )
           tot_charge_ = tot_charge
           !
           ! Make hessian for FCP
           CALL fcp_hessian( helec )
           IF ( capacitance > eps4 ) THEN
              helec = MIN( capacitance, helec )
           END IF
           !
        END IF
        !
        IF ( ANY( if_pos(:,:) == 1 ) .OR. lmovecell .OR. lfcp ) THEN
           !
           filebfgs = TRIM(tmp_dir) // TRIM(prefix) // '.bfgs'
           CALL bfgs( filebfgs, pos, h, relec, etot, grad, fcell, iforceh, &
                      felec, epse, epsf, epsp1, fcp_eps, energy_error, &
                      gradient_error, cell_error, fcp_error, lmovecell, lfcp, &
                      capacitance, helec, step_accepted, conv_ions, &
                      optimizer_failed, istep )
           !
        ELSE
           !
           step_accepted = .FALSE.
           conv_ions     = .TRUE.
           istep         = istep + 1
           !
        END IF
        !
        IF ( lmovecell ) THEN
           ! changes needed only if cell moves
           IF (fix_volume) CALL impose_deviatoric_strain( alat*at, h )
           IF (fix_area)   CALL impose_deviatoric_strain_2d( alat*at, h )
           at = h / alat
           IF(enforce_ibrav) CALL remake_cell( ibrav, alat, at(1,1),at(1,2),at(1,3), new_alat )
           CALL recips( at(1,1),at(1,2),at(1,3), bg(1,1),bg(1,2),bg(1,3) )
           CALL volume( alat, at(1,1),at(1,2),at(1,3), omega )
           !
        ENDIF
        !
        IF ( lfcp ) THEN
           nelec = relec
           tot_charge = SUM(zv(ityp(1:nat))) - nelec
        END IF
        !
        CALL cryst_to_cart( nat, pos,  at, 1 )
        tau    =  RESHAPE( pos,  (/ 3, nat /) )
        !
        IF(enforce_ibrav) CALL output_tau_rescaled(alat/new_alat)
        !
        DEALLOCATE( pos, grad )
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
              ELSEIF ( lmovecell ) THEN
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
              ENDIF
              !
           ELSEIF ( ions_status == 2 ) THEN
              !
              ! ... check with nonzero magnetization succeeded, see above
              !
              IF ( lmovecell ) THEN
                 ions_status = 1
              ELSE
                 ions_status = 0
              ENDIF
              !
           ELSEIF ( ions_status == 1 ) THEN
              !
              ions_status = 0
              !
           ENDIF
           !
           IF ( ions_status < 2 ) THEN
              !
              IF ( ANY( if_pos(:,:) == 1 ) .OR. lmovecell .OR. lfcp ) THEN
                 !
                 CALL terminate_bfgs ( etot, epse, epsf, epsp, fcp_eps, &
                                       lmovecell, lfcp, optimizer_failed )
                 !
              END IF
              !
           END IF
           !
        ELSEIF ( idone == nstep ) THEN
           !
           CALL terminate_bfgs( etot, epse, epsf, epsp, fcp_eps, &
                                lmovecell, lfcp, optimizer_failed )
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
           ENDIF
           !
           IF ( tr2 > 1.D-10 ) THEN
              WRITE( stdout, &
                     '(5X,"new conv_thr",T30,"= ",0PF18.10," Ry",/)' ) tr2
           ELSE
              WRITE( stdout, &
                     '(5X,"new conv_thr",T30,"= ",1PE18.1 ," Ry",/)' ) tr2
           ENDIF
           !
        ENDIF
        !
        CALL output_tau( lmovecell, conv_ions )
        !
        IF ( lfcp ) THEN
           CALL output_fcp( tot_charge_, conv_ions )
        END IF
        !
     ENDIF bfgs_minimization
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
           IF ( ANY( if_pos(:,:) == 1 ) ) THEN
              !
              CALL proj_verlet( conv_ions )
              !
           ELSE
              !
              conv_ions = .TRUE.
              !
           END IF
           !
           ! ... relaxation of FCP
           !
           IF ( lfcp ) THEN
              !
              CALL fcp_relax( conv_ions )
              !
              ! ... finalize FCP
              !
              IF ( conv_ions ) CALL fcp_terminate()
              !
           END IF
           !
           IF ( .NOT. conv_ions .AND. idone >= nstep ) THEN
              !
              WRITE( UNIT = stdout, FMT =  &
                   '(/,5X,"The maximum number of steps has been reached.")' )
              WRITE( UNIT = stdout, &
                   FMT = '(/,5X,"End of molecular dynamics calculation")' )
              !
              ! ... finalize FCP
              !
              IF ( lfcp ) CALL fcp_terminate()
              !
           ENDIF
           !
        ELSEIF ( calc == 'fi' ) THEN
           !
           CALL fire( conv_ions)
           ! 
           IF ( .NOT. conv_ions .AND. idone >= nstep ) THEN
              WRITE( UNIT = stdout, FMT =  &
                   '(/,5X,"The maximum number of steps has been reached.")' )
              WRITE( UNIT = stdout, &
                   FMT = '(/,5X,"End of FIRE minimization")' )
           ENDIF
        ELSEIF ( calc(1:1) == 'l' ) THEN
           !
           ! ... for smart monte carlo method
           !
           IF ( calc(2:2) == 's' ) CALL smart_MC()
           !
           CALL langevin_md()
           !
           ! ... FCP does not support Langevin-dynamics
           !
           IF ( lfcp ) CALL errore('move_ions', &
                         & 'FCP does not support Langevin-dynamics', 1)
           !
           IF ( idone >= nstep ) THEN
              WRITE( UNIT = stdout, FMT =  &
                   '(/,5X,"The maximum number of steps has been reached.")' )
              WRITE( UNIT = stdout, &
                   FMT = '(/,5X,"End of molecular dynamics calculation")' )
              conv_ions = .true.
           ENDIF
           !
        ELSEIF ( calc == 'vd' ) THEN
           !
           IF ( ANY( if_pos(:,:) == 1 ) ) THEN
              !
              IF (tnosep) THEN 
                CALL ions_nosevel(vnhp, xnhp0, xnhpm, RyDt_to_HaDt * dt, nhpcl, nhpdim) 
                ions_nose_energy = ions_nose_nrg(xnhp0,vnhp,qnp,gkbt2nhp,kbt,nhpcl, nhpdim)
              END IF  
              CALL verlet()
              if (idone == 1) nose_cycle = 0
              IF (tnosep) THEN 
                DO 
                  CALL ions_noseupd(xnhpp, xnhp0, xnhpm, RyDt_to_HaDt * dt, qnp, ekin2nhp, gkbt2nhp, vnhp, kbt, &
                               nhpcl, nhpdim, nhpbeg, nhpend)
                  CALL ions_nose_shiftvar(xnhpp, xnhp0, xnhpm)
                  nose_cycle = nose_cycle + 1
                  IF (nose_cycle .ge. 2) EXIT 
                END DO 
                
              END IF         
              !
           END IF
           !
           ! ... the dynamics of FCP
           !
           IF ( lfcp ) CALL fcp_verlet()
           !
           IF ( idone >= nstep) THEN
              !
              CALL terminate_verlet()
              !
              ! ... finalize FCP
              !
              IF ( lfcp ) CALL fcp_terminate()
              !
              conv_ions = .true.
              !
           ENDIF
        ELSE IF (calc .eq. 'wd' .AND. ANY(if_pos(:,:) == 1) ) THEN
            CALL velocity_verlet() 
            IF (idone .GE. nstep) THEN 
               CALL terminate_verlet() 
               conv_ions = .true.  
            END IF
        ELSE
           !
           ! ... variable cell shape md
           !
           CALL vcsmd( conv_ions )
           !
           ! ... FCP does not support cell shape md
           !
           IF ( lfcp ) CALL errore('move_ions', &
                         & 'FCP does not support cell shape MD', 1)
           !
           ! ... after nstep, set conv_ions to T for MD, to F for damped MD
           !
           IF ( .NOT.conv_ions .AND. idone >= nstep ) THEN
              WRITE( UNIT = stdout, FMT = '(/,5X,"Maximum number of ", &
             &    "iterations reached, stopping")' )
              conv_ions = ( calc(2:2) == 'd' )
           ENDIF
           !
        ENDIF
        !
        IF ( conv_ions ) ions_status  = 0
        !
     ENDIF
     !
     ! ... before leaving check that the new positions still transform
     ! ... according to the symmetry of the system.
     ! ... FIXME: should be done in all cases, not just for vc-md
     ! ... FIXME 2: why not impose symmetry instead of just checking it?
     !
     CALL checkallsym( nat, tau, ityp)

     ! write trajectory output files

     if (mod(istep, iprint)==0) then
        iunit = printout_base_open('.pos')
        call printout_pos(iunit, tau*alat, nat, tps=elapsed_time, nfi=istep)
        call printout_base_close(iunit)
        iunit = printout_base_open('.cel')
        call printout_cell(iunit,at*alat,istep,elapsed_time)
        call printout_base_close(iunit)
        iunit = printout_base_open('.for')
        call printout_pos(iunit, force*alat, nat, tps=elapsed_time, nfi=istep)
        call printout_base_close(iunit)
        iunit = printout_base_open('.vel')
        call printout_pos(iunit, vel*alat, nat, tps=elapsed_time, nfi=istep)
        call printout_base_close(iunit)
     endif


     !
  ENDIF
  !
  !
  CALL mp_bcast( ions_status, ionode_id, intra_image_comm )
  CALL mp_bcast( optimizer_failed, ionode_id, intra_image_comm )
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
  ENDIF
  !
  IF ( lfcp ) THEN
     CALL mp_bcast(nelec,      ionode_id, intra_image_comm)
     CALL mp_bcast(tot_charge, ionode_id, intra_image_comm)
  END IF
  !
  !
  ! ... update convergence threshold of 3D-RISM
  !
  IF ( lrism ) THEN
     IF ( tr2 < starting_scf_threshold ) THEN
       CALL rism_new_conv_thr()
     END IF
  END IF
  !
  RETURN
  !
END SUBROUTINE move_ions
