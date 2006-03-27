!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE cprmain( tau, fion_out, etot_out )
  !----------------------------------------------------------------------------
  !
  USE kinds,                    ONLY : DP
  USE constants,                ONLY : bohr_radius_angs, amu_au
  USE control_flags,            ONLY : iprint, isave, thdyn, tpre, tbuff,      &
                                       iprsta, tfor, tvlocw,      &
                                       taurdr, tprnfor, tsdc, lconstrain, lwf, &
                                       lneb, lcoarsegrained, ndr, ndw, nomore, &
                                       tsde, tortho, tnosee, tnosep, trane,    &
                                       tranp, tsdp, tcp, tcap, ampre, amprp,   &
                                       tnoseh, tolp, ortho_eps, ortho_max,     &
                                       printwfc
  USE core,                     ONLY : nlcc_any, rhoc
  USE uspp_param,               ONLY : nhm, nh
  USE cvan,                     ONLY : nvb, ish
  USE uspp,                     ONLY : nkb, vkb, becsum, deeq
  USE energies,                 ONLY : eht, epseu, exc, etot, eself, enl, &
                                       ekin, atot, entropy, egrand, enthal, &
                                       ekincm, print_energies
  USE electrons_base,           ONLY : nbspx, nbsp, ispin, f, nspin
  USE electrons_base,           ONLY : nel, iupdwn, nupdwn, nudx, nelt
  USE efield_module,            ONLY : efield, epol, tefield, allocate_efield, &
                                       efield_update, ipolp, qmat, gqq, evalue,&
                                       berry_energy, pberryel, pberryion,      &
                                       efield2, epol2, tefield2,               &
                                       allocate_efield2, efield_update2,       &
                                       ipolp2, qmat2, gqq2, evalue2,           &
                                       berry_energy2, pberryel2, pberryion2
  USE ensemble_dft,             ONLY : tens, tgrand, ninner, ismear, etemp,   &
                                       ef, tdynz, tdynf, zmass, fmass, fricz, &
                                       fricf, allocate_ensemble_dft,          &
                                       id_matrix_init, z0, c0diag, becdiag,   &
                                       bec0, v0s, vhxcs, becdrdiag, gibbsfe
  USE cg_module,                ONLY : tcg, maxiter, etresh, passop, &
                                       allocate_cg, cg_update, &
                                       itercg, c0old
  USE gvecp,                    ONLY : ngm
  USE gvecs,                    ONLY : ngs
  USE gvecb,                    ONLY : ngb
  USE gvecw,                    ONLY : ngw
  USE reciprocal_vectors,       ONLY : gstart, mill_l
  USE ions_base,                ONLY : na, nat, pmass, nax, nsp, rcmax
  USE ions_base,                ONLY : ind_srt, ions_cofmass, ions_kinene, &
                                       ions_temp, ions_thermal_stress, if_pos
  USE ions_base,                ONLY : ions_vrescal, fricp, greasp, &
                                       iforce, ndfrz, ions_shiftvar, ityp, &
                                       atm, ind_bck, cdm, cdms
  USE cell_base,                ONLY : a1, a2, a3, b1, b2, b3, ainv, frich, &
                                       greash, tpiba2, omega, alat, ibrav,  &
                                       celldm, h, hold, hnew, velh, deth,   &
                                       wmass, press, iforceh, cell_force,   &
                                       thdiag
  USE grid_dimensions,          ONLY : nnrx, nr1, nr2, nr3
  USE smooth_grid_dimensions,   ONLY : nnrsx, nr1s, nr2s, nr3s
  USE smallbox_grid_dimensions, ONLY : nr1b, nr2b, nr3b
  USE local_pseudo,             ONLY : allocate_local_pseudo
  USE io_global,                ONLY : io_global_start, stdout, ionode
  USE dener,                    ONLY : detot
  USE derho,                    ONLY : drhor, drhog
  USE cdvan,                    ONLY : dbec, drhovan
  USE stre,                     ONLY : stress
  USE gvecw,                    ONLY : ggp
  USE parameters,               ONLY : nsx
  USE constants,                ONLY : pi, factem, au_gpa, au_ps, gpa_au
  USE io_files,                 ONLY : psfile, pseudo_dir
  USE wave_base,                ONLY : wave_steepest, wave_verlet
  USE wave_base,                ONLY : wave_speed2, frice, grease
  USE control_flags,            ONLY : conv_elec, tconvthrs
  USE check_stop,               ONLY : check_stop_now
  USE efcalc,                   ONLY : clear_nbeg, ef_force
  USE ions_base,                ONLY : zv, ions_vel
  USE cp_electronic_mass,       ONLY : emass, emass_cutoff, emass_precond
  USE ions_positions,           ONLY : tau0, taum, taup, taus, tausm, tausp, &
                                       vels, velsm, velsp, ions_hmove,       &
                                       ions_move, fion, fionm
  USE ions_nose,                ONLY : gkbt, kbt, qnp, ndega, nhpcl, nhpdim, &
                                       nhpend, vnhp, xnhp0, xnhpm, xnhpp,    &
                                       atm2nhp, ions_nosevel, ions_noseupd,  &
                                       tempw, ions_nose_nrg, gkbt2nhp,       &
                                       ekin2nhp, anum2nhp
  USE electrons_nose,           ONLY : qne, ekincw, xnhe0, xnhep, xnhem,  &
                                       vnhe, electrons_nose_nrg,    &
                                       electrons_nose_shiftvar,           &
                                       electrons_nosevel, electrons_noseupd
  USE from_scratch_module,      ONLY : from_scratch
  USE from_restart_module,      ONLY : from_restart
  USE wavefunctions_module,     ONLY : c0, cm, phi => cp
  USE wannier_module,           ONLY : allocate_wannier
  USE print_out_module,         ONLY : printout_new
  USE printout_base,            ONLY : printout_base_open, &
                                       printout_base_close, &
                                       printout_pos, printout_cell, &
                                       printout_stress
  USE cell_nose,                ONLY : xnhh0, xnhhm, xnhhp, vnhh, temph, &
                                       qnh, cell_nosevel, cell_noseupd,  &
                                       cell_nose_nrg, cell_nose_shiftvar
  USE cell_base,                ONLY : cell_kinene, cell_gamma, &
                                       cell_move, cell_hmove
  USE gvecw,                    ONLY : ecutw
  USE gvecp,                    ONLY : ecutp
  USE time_step,                ONLY : delt, tps, dt2, dt2by2, twodelt
  USE electrons_module,         ONLY : cp_eigs
  USE print_out_module,         ONLY : cp_print_rho
  USE cp_main_variables,        ONLY : allocate_mainvar, &
                                       acc, bec, lambda, lambdam, lambdap, &
                                       ema0bg, sfac, eigr, ei1, ei2, ei3,  &
                                       irb, becdr, taub, eigrb, rhog, rhos, &
                                       rhor, rhopr, bephi, becp, nfi
  USE autopilot,                ONLY : event_step, event_index, &
                                       max_event_step, restart_p
  USE metadyn_vars,             ONLY : dfe_acc, etot_av
  USE cell_base,                ONLY : s_to_r, r_to_s
  USE phase_factors_module,     ONLY : strucf
  USE cpr_subroutines,          ONLY : print_lambda, print_atomic_var, &
                                       ions_cofmsub, elec_fakekine
  USE wannier_subroutines,      ONLY : wannier_startup, wf_closing_options, &
                                       ef_enthalpy
  USE restart_file,             ONLY : readfile, writefile
  USE constraints_module,       ONLY : check_constraint, lagrange, &
                                       remove_constr_force
  USE metadyn_base,             ONLY : set_target
  USE autopilot,                ONLY : pilot
  USE ions_nose,                ONLY : ions_nose_allocate, ions_nose_shiftvar
  USE orthogonalize,            ONLY : ortho
  USE orthogonalize_base,       ONLY : updatc
  USE control_flags,            ONLY : force_pairing
  !
  IMPLICIT NONE
  !
  ! ... input/output variables
  !
  REAL(DP), INTENT(INOUT) :: tau(3,nat)
  REAL(DP), INTENT(OUT)   :: fion_out(3,nat)
  REAL(DP), INTENT(OUT)   :: etot_out
  !
  ! ... control variables
  !
  LOGICAL :: tfirst, tlast, tstop, tconv
  LOGICAL :: ttprint    
    !  logical variable used to control printout
  !
  ! ... forces on ions
  !
  REAL(DP) :: maxfion
  !
  ! ... work variables
  !
  REAL(DP) :: tempp, savee, saveh, savep, epot, epre, &
              enow, econs, econt, fccc, ccc, bigr, dt2bye
  REAL(DP) :: ekinc0, ekinp, ekinpr, ekinc
  REAL(DP) :: temps(nsx)
  REAL(DP) :: ekinh, temphc, randy
  REAL(DP) :: delta_etot
  REAL(DP) :: ftmp, enb, enbi
  INTEGER  :: is, nacc, ia, j, iter, i, isa, ipos
  INTEGER   :: k, ii, l, m, iss
  !
  REAL(DP) :: hgamma(3,3), temphh(3,3)
  REAL(DP) :: fcell(3,3)
  !
  REAL(DP) :: stress_gpa(3,3), thstress(3,3)
  !
  REAL(DP), ALLOCATABLE :: tauw(:,:)  
    ! temporary array used to printout positions
  CHARACTER(LEN=3) :: labelw( nat )
  ! for force_pairing
  INTEGER   :: n_spin_start 
  !
  dt2bye   = dt2 / emass
  etot_out = 0.D0
  enow     = 1.D9
  !
  tfirst = .TRUE.
  tlast  = .FALSE.
  nacc   = 5
  !
  n_spin_start = nspin
  IF( force_pairing ) n_spin_start = 1
  !
  ! ... Check for restart_p from Autopilot Feature Suite
  !
  IF ( restart_p ) THEN
     !
     ! ... do not add past nfi
     !
     nomore = nomore
     !
  END IF
  !
  !======================================================================
  !
  !           basic loop for molecular dynamics starts here
  !
  !======================================================================
  !
  main_loop: DO
     !
     CALL start_clock( 'total_time' )
     !
     nfi     = nfi + 1
     tlast   = ( nfi == nomore )
     ttprint = ( MOD( nfi, iprint ) == 0 ).or.tlast
     !
     IF ( ionode .AND. ttprint ) &
        WRITE( stdout, '(/," * Physical Quantities at step:",I6)' ) nfi
     !
     IF ( tsde ) THEN
        fccc = 1.D0 
     ELSE
        fccc = 1.D0 / ( 1.D0 + frice )
     END IF
     !
     ! ... calculation of velocity of nose-hoover variables
     !
     IF ( tnosep ) CALL ions_nosevel( vnhp, xnhp0, xnhpm, delt, nhpcl, nhpdim )
     !
     IF ( tnosee ) THEN
        !
        CALL electrons_nosevel( vnhe, xnhe0, xnhem, delt )
        !
        fccc = 1.D0 / ( 1.D0 + 0.5D0 * delt * vnhe )
        !
     END IF
     !
     IF ( tnoseh ) THEN
        !
        CALL cell_nosevel( vnhh, xnhh0, xnhhm, delt )
        !
        velh(:,:) = 2.D0 * ( h(:,:) - hold(:,:) ) / delt - velh(:,:)
        !
     END IF
     ! 
     IF ( tfor .OR. thdyn .OR. tfirst ) THEN
        !
        CALL initbox( tau0, taub, irb )
        CALL phbox( taub, eigrb )
        !
     END IF
     !
     IF ( tfor .OR. thdyn ) CALL phfac( tau0, ei1, ei2, ei3, eigr ) 
     !
     ! ... strucf calculates the structure factor sfac
     !
     CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
     !
     IF ( thdyn ) CALL formf( tfirst, eself )
     !
     ! ... why this call ??? from Paolo Umari
     !
     IF ( tefield.or.tefield2 ) CALL calbec( 1, nsp, eigr, c0, bec ) ! ATTENZIONE  
     !
     ! Autopilot (Dynamic Rules) Implimentation    
     !
     call pilot(nfi)
    
     !
     IF ( ( tfor .OR. tfirst ) .AND. tefield ) CALL efield_update( eigr )
     IF ( ( tfor .OR. tfirst ) .AND. tefield2 ) CALL efield_update2( eigr )

     !
     !=======================================================================
     !
     !    electronic degrees of freedom are updated here
     !
     !=======================================================================
     !
     IF( force_pairing ) THEN
          c0(:,iupdwn(2):nbsp,1,1)       =     c0(:,1:nupdwn(2),1,1)
          cm(:,iupdwn(2):nbsp,1,1)       =     cm(:,1:nupdwn(2),1,1)
         phi(:,iupdwn(2):nbsp,1,1)       =    phi(:,1:nupdwn(2),1,1)
      lambda(1:nupdwn(2),1:nupdwn(2), 2) = lambda(1:nupdwn(2),1:nupdwn(2), 1)
      lambda(nudx, nudx, 2) = 0.d0
     ENDIF
     !
     !
     CALL move_electrons( nfi, tfirst, tlast, b1, b2, b3, fion, &
                          enthal, enb, enbi, fccc, ccc, dt2bye )
     !
     IF ( tpre ) THEN
        !
        CALL nlfh( bec, dbec, lambda )
        !
        CALL ions_thermal_stress( stress, pmass, omega, h, vels, nsp, na )
        !
     END IF
     !
     !=======================================================================
     !
     !              verlet algorithm
     !
     !     loop which updates cell parameters and ionic degrees of freedom
     !     hnew=h(t+dt) is obtained from hold=h(t-dt) and h=h(t)
     !     tausp=pos(t+dt) from tausm=pos(t-dt) taus=pos(t) h=h(t)
     !
     !           guessed displacement of ions
     !=======================================================================
     !
     hgamma(:,:) = 0.D0
     !
     IF ( thdyn ) THEN
        !
        CALL cell_force( fcell, ainv, stress, omega, press, wmass )
        !
        CALL cell_move( hnew, h, hold, delt, iforceh, &
                        fcell, frich, tnoseh, vnhh, velh, tsdc )
        !
        velh(:,:) = ( hnew(:,:) - hold(:,:) ) / twodelt
        !
        CALL cell_gamma( hgamma, ainv, h, velh )
        !
     END IF
     !
     !======================================================================
     !
     IF ( tfor ) THEN
        !
        IF ( lwf ) CALL ef_force( fion, na, nsp, zv )
        !
        IF ( lconstrain ) THEN
           !
           IF ( lcoarsegrained ) CALL set_target()
           !
           ! ... we first remove the component of the force along the 
           ! ... constrain gradient (this constitutes the initial guess 
           ! ... for the lagrange multiplier)
           !
           CALL remove_constr_force( nat, tau0, &
                                     iforce, ityp, 1.D0, fion(:3,:nat) )
           !
        END IF
        !
        CALL ions_move( tausp, taus, tausm, iforce, pmass, fion, ainv, &
                        delt, na, nsp, fricp, hgamma, vels, tsdp, tnosep, &
                        fionm, vnhp, velsp, velsm, nhpcl, nhpdim, atm2nhp )
        !
        IF ( lconstrain ) THEN
           !
           ! ... constraints are imposed here
           !
           CALL s_to_r( tausp, taup, na, nsp, hnew )
           !
           CALL check_constraint( nat, taup, tau0, fion(:3,:nat), &
                                  iforce, ityp, 1.D0, delt, amu_au )
           !
           CALL r_to_s( taup, tausp, na, nsp, ainv )
           !
           ! ... average value of the lagrange multipliers
           !
           IF ( lcoarsegrained ) THEN
              !
              etot_av = etot_av + etot
              !
              dfe_acc(:) = dfe_acc(:) - lagrange(:)
              !
           END IF
           !
        END IF
        !
        CALL ions_cofmass( tausp, pmass, na, nsp, cdm )
        !
        IF ( ndfrz == 0 ) &
           CALL ions_cofmsub( tausp, iforce, na, nsp, cdm, cdms )
        !
        CALL s_to_r( tausp, taup, na, nsp, hnew )
        !
     END IF
     !     
     !--------------------------------------------------------------------------
     !              initialization with guessed positions of ions
     !--------------------------------------------------------------------------
     !
     ! ... if thdyn=true g vectors and pseudopotentials are recalculated for 
     ! ... the new cell parameters
     !
     IF ( tfor .OR. thdyn ) THEN
        !
        IF ( thdyn ) THEN
           !
           hold = h
           h    = hnew
           !
           CALL newinit( h )
           !
           CALL newnlinit()
           !
        ELSE
           !
           hold = h
           !
        END IF
        !
        ! ... phfac calculates eigr
        !
        CALL phfac( taup, ei1, ei2, ei3, eigr )
        !
        ! ... prefor calculates vkb
        !
        CALL prefor( eigr, vkb )
        !
     END IF
     !
     !--------------------------------------------------------------------------
     !                    imposing the orthogonality
     !--------------------------------------------------------------------------
     !
     IF ( .NOT. tcg ) THEN
        !
        IF ( tortho ) THEN
           !
           CALL ortho( eigr, cm(:,:,1,1), phi(:,:,1,1), &
                       lambda, bigr, iter, ccc, bephi, becp )
           !
        ELSE
           !
           CALL gram( vkb, bec, nkb, cm, ngw, nbsp )
           !
           IF ( iprsta > 4 ) CALL dotcsc( eigr, cm )
           !
        END IF
        !
     END IF
     !
     !--------------------------------------------------------------------------
     !                   correction to displacement of ions
     !--------------------------------------------------------------------------
     !
     IF ( .NOT. tcg ) THEN
        !
        IF ( iprsta >= 3 ) CALL print_lambda( lambda, nbsp, 9, 1.D0 )
        !
        IF ( tortho ) THEN
           DO iss = 1, n_spin_start
              CALL updatc( ccc, nbsp, lambda(:,:,iss), SIZE(lambda,1), phi, SIZE(phi,1), &
                        bephi, SIZE(bephi,1), becp, bec, cm, nupdwn(iss), iupdwn(iss) )
           END DO
        END IF
        !
        IF( force_pairing ) THEN
              c0(:,iupdwn(2):nbsp,1,1)       =     c0(:,1:nupdwn(2),1,1)
              cm(:,iupdwn(2):nbsp,1,1)       =     cm(:,1:nupdwn(2),1,1)
             phi(:,iupdwn(2):nbsp,1,1)       =    phi(:,1:nupdwn(2),1,1)
          lambda(1:nupdwn(2),1:nupdwn(2), 2) = lambda(1:nupdwn(2),1:nupdwn(2), 1)
          lambda(nudx, nudx, 2) = 0.d0
        ENDIF
        !
        CALL calbec( nvb+1, nsp, eigr, cm, bec )
        !
        IF ( tpre ) &
           CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, cm, dbec, .TRUE. )
        !
        IF ( iprsta >= 3 ) CALL dotcsc( eigr, cm )
        !
     END IF
     !
     !--------------------------------------------------------------------------
     !                  temperature monitored and controlled
     !--------------------------------------------------------------------------
     !
     ekinp  = 0.D0
     ekinpr = 0.D0
     tempp  = 0.D0
     temps  = 0.D0
     ekinc0 = 0.0d0
     ekinc = 0.0d0
     !
     !
     ! ... ionic kinetic energy 
     !
     IF ( tfor ) THEN
        !
        CALL ions_vel( vels, tausp, tausm, na, nsp, delt )
        !
        CALL ions_kinene( ekinp, vels, na, nsp, hold, pmass )
        !
     END IF
     !
     ! ... ionic temperature
     !
     IF ( tfor ) &
        CALL ions_temp( tempp, temps, ekinpr, vels, na, nsp, &
                        hold, pmass, ndega, nhpdim, atm2nhp, ekin2nhp )
     !
     ! ... fake electronic kinetic energy
     !
     IF ( .NOT. tcg ) THEN
        !
        CALL elec_fakekine( ekinc0, ema0bg, emass, c0, cm, ngw, nbsp, delt )
        !
        ekinc = ekinc0
        !
     END IF
     
     !
     ! ... fake cell-parameters kinetic energy
     !
     ekinh = 0.D0
     !
     IF ( thdyn ) CALL cell_kinene( ekinh, temphh, velh )
     !
     IF ( COUNT( iforceh == 1 ) > 0 ) THEN
        !
        temphc = 2.D0 * factem * ekinh / DBLE( COUNT( iforceh == 1 ) )
        !
     ELSE
        !
        temphc = 0.D0
        !
     END IF
     !
     ! ... udating nose-hoover friction variables
     !
     IF ( tnosep ) CALL ions_noseupd( xnhpp, xnhp0, xnhpm, delt, qnp, &
                                      ekin2nhp, gkbt2nhp, vnhp, kbt,  &
                                      nhpcl, nhpdim, nhpend )
     !
     IF ( tnosee ) CALL electrons_noseupd( xnhep, xnhe0, xnhem, &
                                           delt, qne, ekinc, ekincw, vnhe )
     !
     IF ( tnoseh ) CALL cell_noseupd( xnhhp, xnhh0, xnhhm, &
                                      delt, qnh, temphh, temph, vnhh )
     !
     ! ... warning:  thdyn and tcp/tcap are not compatible yet!!!
     !
     IF ( tcp .OR. tcap .AND. tfor .AND. .NOT.thdyn ) THEN
        !
        IF ( tempp > (tempw+tolp) .OR. &
             tempp < (tempw-tolp) .AND. tempp /= 0.D0 ) THEN
           !
           CALL  ions_vrescal( tcap, tempw, tempp, taup, &
                               tau0, taum, na, nsp, fion, iforce, pmass, delt )
           !
        END IF
        !
     END IF
     !
     IF ( MOD( nfi, iprint ) == 0 .OR. tlast ) THEN
        !
        IF( force_pairing )  THEN
           lambda(:, :, 2) =  lambda(:, :, 1)
           lambdap(:, :, 2) = lambdap(:, :, 1)
           lambda(nudx, nudx, 2) = 0.d0  
           lambdap(nudx, nudx, 2) = 0.d0
           WRITE( stdout, '("Occupations in CPR:")' )
           WRITE( stdout, '(10F9.6)' ) ( f(i), i = 1, nbspx )  
        END IF
        !
        CALL cp_eigs( nfi, lambdap, lambda )
        !
        !
     END IF
     !
     IF ( lwf ) CALL ef_enthalpy( enthal, tau0 )
     !
     IF ( tens ) THEN
        !
        IF ( MOD( nfi, iprint ) == 0 .OR. tlast ) THEN
           !
           WRITE( stdout, '("Occupations  :")' )
           WRITE( stdout, '(10F9.6)' ) ( f(i), i = 1, nbsp )
           !
        END IF
        !
     END IF
     !
     epot = eht + epseu + exc
     !
     acc(1) = acc(1) + ekinc
     acc(2) = acc(2) + ekin
     acc(3) = acc(3) + epot
     acc(4) = acc(4) + etot
     acc(5) = acc(5) + tempp
     !
     IF ( .NOT. tcg ) THEN
        !
        econs = ekinp + ekinh + enthal
        econt = econs + ekinc
        !
     ELSE
        !
        IF ( .NOT. tens ) THEN
           !
           econs = ekinp + etot
           atot  = etot
           econt = econs
           !
        ELSE
           !
           gibbsfe = atot
           econs   = ekinp + atot
           econt   = econs
           !
        END IF
        !
     END IF
     !
     ! ... add energies of thermostats
     !
     IF ( tnosep ) &
        econt = econt + ions_nose_nrg( xnhp0, vnhp, qnp, &
                                       gkbt2nhp, kbt, nhpcl, nhpdim )
     IF ( tnosee ) &
        econt = econt + electrons_nose_nrg( xnhe0, vnhe, qne, ekincw )
     IF ( tnoseh ) &
        econt = econt + cell_nose_nrg( qnh, xnhh0, vnhh, temph, iforceh )
     !
     tps = tps + delt * au_ps
     !
     CALL printout_new( nfi, tfirst, ttprint, ttprint, tps, hold, stress, &
                        tau0, vels, fion, ekinc, temphc, tempp, temps, etot, &
                        enthal, econs, econt, vnhh, xnhh0, vnhp, xnhp0, atot )
     !
     IF( tfor ) THEN
        !
        ! ... new variables for next step
        !
        CALL ions_shiftvar( taup,  tau0, taum  )   !  real positions
        CALL ions_shiftvar( tausp, taus, tausm )   !  scaled positions         
        CALL ions_shiftvar( velsp, vels, velsm )   !  scaled velocities
        !
        IF ( tnosep ) CALL ions_nose_shiftvar( xnhpp, xnhp0, xnhpm )
        IF ( tnosee ) CALL electrons_nose_shiftvar( xnhep, xnhe0, xnhem )
        IF ( tnoseh ) CALL cell_nose_shiftvar( xnhhp, xnhh0, xnhhm )
        !
     END IF
     !
     IF ( thdyn ) CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emass_cutoff )
     !
     ekincm = ekinc0
     !  
     ! ... cm=c(t+dt) c0=c(t)
     !
     IF( .NOT. tcg ) THEN
        !
        CALL dswap( 2*ngw*nbsp, c0, 1, cm, 1 )
        !
     ELSE
        !
        CALL cg_update( tfirst, nfi, c0 )
        !
        IF ( tfor .AND. .NOT. tens .AND. &
             ( ( MOD( nfi, isave ) == 0 ) .OR. tlast ) ) THEN
           !
           ! ... in this case optimize c0 and lambda for smooth
           ! ... restart with CP
           !
           CALL initbox( tau0, taub, irb )
           CALL phbox( taub, eigrb )
           CALL phfac( tau0, ei1, ei2, ei3, eigr )
           CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
           !
           IF ( thdyn )    CALL formf( tfirst, eself )
           IF ( tefield )  CALL efield_update( eigr )
           IF ( tefield2 ) CALL efield_update2( eigr )
           !
           lambdam = lambda
           !
           CALL move_electrons( nfi, tfirst, tlast, b1, b2, b3, &
                                fion, enthal, enb, enbi, fccc, ccc, dt2bye )
           !
        END IF
        !
     END IF
     !
     ! ... now:  cm=c(t) c0=c(t+dt)
     !
     tfirst = .FALSE.
     !
     ! ... write on file ndw each isave
     !
     IF ( ( MOD( nfi, isave ) == 0 ) .AND. ( nfi < nomore ) ) THEN
        !
        IF ( tcg ) THEN
          !
          CALL writefile( ndw, h, hold ,nfi, c0(:,:,1,1), c0old, taus, tausm,  &
                          vels, velsm, acc, lambda, lambdam, xnhe0, xnhem,     &
                          vnhe, xnhp0, xnhpm, vnhp, nhpcl,nhpdim,ekincm, xnhh0,&
                          xnhhm, vnhh, velh, ecutp, ecutw, delt, pmass, ibrav, &
                          celldm, fion, tps, z0, f, rhopr )
           !
        ELSE
           !
           CALL writefile( ndw, h, hold, nfi, c0(:,:,1,1), cm(:,:,1,1), taus,  &
                           tausm, vels, velsm, acc,  lambda, lambdam, xnhe0,   &
                           xnhem, vnhe, xnhp0, xnhpm, vnhp,nhpcl,nhpdim,ekincm,&
                           xnhh0, xnhhm, vnhh, velh, ecutp, ecutw, delt, pmass,&
                           ibrav, celldm, fion, tps, z0, f, rhopr )
           !
        END IF
        !
     END IF
     !
     epre = enow
     enow = etot
     !
     frice = frice * grease
     fricp = fricp * greasp
     frich = frich * greash
     !
     !======================================================================
     !
     CALL stop_clock( 'total_time' )
     !
     delta_etot = ABS( epre - enow )
     !
     tstop = check_stop_now()
     tconv = .FALSE.
     !
     IF ( tconvthrs%active ) THEN
        !
        ! ... electrons
        !
        tconv = ( ekinc < tconvthrs%ekin .AND. delta_etot < tconvthrs%derho )
        !
        IF ( tfor ) THEN
           !
           ! ... ions
           !
           maxfion = MAXVAL( ABS( fion(:,1:nat) ) )
           !
           tconv = tconv .AND. ( maxfion < tconvthrs%force )
           !
        END IF
        !
     END IF
     !
     ! ... in the case cp-wf the check on convergence is done starting
     ! ... from the second step 
     !
     IF ( lwf .AND. tfirst ) tconv = .FALSE.
     !
     IF ( tconv ) THEN
        !
        IF ( ionode ) THEN
           !
           WRITE( stdout, &
                & "(/,3X,'MAIN:',10X,'EKINC   (thr)', &
                & 10X,'DETOT   (thr)',7X,'MAXFORCE   (thr)')" )
           WRITE( stdout, "(3X,'MAIN: ',3(D14.6,1X,D8.1))" ) &
               ekinc, tconvthrs%ekin, delta_etot,                  &
               tconvthrs%derho, 0.D0, tconvthrs%force
           WRITE( stdout, &
                  "(3X,'MAIN: convergence achieved for system relaxation')" )
           !
        END IF
        !
     END IF
     !
     tstop = tstop .OR. tconv
     !
     IF ( lwf ) &
        CALL wf_closing_options( nfi, c0, cm, bec, becdr, eigr, eigrb, taub, &
                                 irb, ibrav, b1, b2, b3, taus, tausm, vels,  &
                                 velsm, acc, lambda, lambdam, xnhe0, xnhem,  &
                                 vnhe, xnhp0, xnhpm, vnhp, nhpcl, nhpdim,    &
                                 ekincm, xnhh0, xnhhm, vnhh, velh, ecutp,    &
                                 ecutw, delt, celldm, fion, tps, z0, f, rhopr )
     !
     IF ( ( nfi >= nomore ) .OR. tstop ) EXIT main_loop
     !
  END DO main_loop
  !
  !===================== end of main loop of molecular dynamics ===============
  ! 
  ! ... Here copy relevant physical quantities into the output arrays/variables
  !
  etot_out = etot
  !
  isa = 0
  !
  DO is = 1, nsp
     !
     DO ia = 1, na(is)
        !
        isa = isa + 1
        !
        ipos = ind_srt( isa )
        !
        tau(:,ipos) = tau0(:,isa)
        !
        fion_out(:,ipos) = fion(:,isa)
        !
     END DO
     !
  END DO
  !
  IF ( lneb ) fion_out(:,1:nat) = fion(:,1:nat) * DBLE( if_pos(:,1:nat) )
  !
  conv_elec = .TRUE.
  !
  acc = acc / DBLE( nfi )
  !
  IF ( tcg ) THEN
     !
     CALL writefile( ndw, h, hold, nfi, c0(:,:,1,1), c0old, taus, tausm, vels, &
                     velsm, acc, lambda, lambdam, xnhe0, xnhem, vnhe, xnhp0,   &
                     xnhpm, vnhp, nhpcl,nhpdim,ekincm, xnhh0, xnhhm, vnhh,     &
                     velh, ecutp, ecutw, delt, pmass, ibrav, celldm, fion,     &
                     tps, z0, f, rhopr )
     !
  ELSE
     !
     CALL writefile( ndw, h, hold, nfi, c0(:,:,1,1), cm(:,:,1,1), taus, tausm, &
                     vels, velsm, acc, lambda, lambdam, xnhe0, xnhem, vnhe,    &
                     xnhp0, xnhpm, vnhp, nhpcl,nhpdim,ekincm, xnhh0, xnhhm,    &
                     vnhh, velh, ecutp, ecutw, delt, pmass, ibrav, celldm,     &
                     fion, tps, z0, f, rhopr )
     !
  END IF
  !
  IF( iprsta > 1 ) CALL print_lambda( lambda, nbsp, nbsp, 1.D0 )
  !
  RETURN
  !
END SUBROUTINE cprmain
!
!----------------------------------------------------------------------------
SUBROUTINE terminate_run()
  !----------------------------------------------------------------------------
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout, ionode
  USE cp_main_variables, ONLY : acc
  !
  IMPLICIT NONE
  !
  INTEGER  :: i, nacc = 5
  !
  ! ...  print statistics
  !
  IF ( ionode ) THEN
     !
     WRITE( stdout, 1949 )
     WRITE( stdout, 1950 ) ( acc(i), i = 1, nacc )
     !
  END IF
  !
1949 FORMAT( //'              averaged quantities :',/,9X, &
           & 'ekinc',10X,'ekin',10X,'epot',10X,'etot',5X,'tempp' )
1950 FORMAT( 4F14.5,F10.1 )
  !
  CALL print_clock( 'initialize' )
  CALL print_clock( 'total_time' )
  CALL print_clock( 'formf' )
  CALL print_clock( 'rhoofr' )
  CALL print_clock( 'vofrho' )
  CALL print_clock( 'dforce' )
  CALL print_clock( 'calphi' )
  CALL print_clock( 'ortho' )
  CALL print_clock( 'updatc' )
  CALL print_clock( 'gram' )
  CALL print_clock( 'newd' )
  CALL print_clock( 'calbec' )
  CALL print_clock( 'prefor' )
  CALL print_clock( 'strucf' )
  CALL print_clock( 'nlfl' )
  CALL print_clock( 'nlfq' )
  CALL print_clock( 'set_cc' )
  CALL print_clock( 'rhov' )
  CALL print_clock( 'nlsm1' )
  CALL print_clock( 'nlsm2' )
  CALL print_clock( 'forcecc' )
  CALL print_clock( 'fft' )
  CALL print_clock( 'ffts' )
  CALL print_clock( 'fftw' )
  CALL print_clock( 'fftb' )
  CALL print_clock( 'rsg' )
  !
1974 FORMAT( 1X,2I5,3F10.4,2X,3F10.4 )
1975 FORMAT( /1X,'Scaled coordinates '/1X,'species',' atom #' )
1976 FORMAT( 1X,2I5,3F10.4 )
  !
  IF ( ionode ) &
     WRITE( stdout, '(5X,//,24("=")," end cp ",24("="),//)' ) 
  !
  CALL memory()
  !
END SUBROUTINE terminate_run
