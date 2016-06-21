!
! Copyright (C) 2002-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE cprmain( tau_out, fion_out, etot_out )
  !----------------------------------------------------------------------------
  !
  USE kinds,                    ONLY : DP
  USE constants,                ONLY : bohr_radius_angs, amu_au, au_gpa
  USE control_flags,            ONLY : iprint, isave, thdyn, tpre, iverbosity, &
                                       tfor, remove_rigid_rot, taurdr, llondon,&
                                       tprnfor, tsdc, lconstrain, lwf,         &
                                       ndr, ndw, nomore, tsde, textfor,        &
                                       tortho, tnosee, tnosep, trane, tranp,   &
                                       tsdp, tcp, tcap, ampre, amprp, tnoseh,  &
                                       tolp, ortho_eps, ortho_max
  USE core,                     ONLY : rhoc
  USE uspp_param,               ONLY : nhm, nh, nvb, ish
  USE uspp,                     ONLY : nkb, vkb, becsum, deeq, okvan, nlcc_any
  USE energies,                 ONLY : eht, epseu, exc, etot, eself, enl, &
                                       ekin, atot, entropy, egrand, enthal, &
                                       ekincm, print_energies
  USE electrons_base,           ONLY : nbspx, nbsp, ispin, f, nspin, nbsp_bgrp
  USE electrons_base,           ONLY : nel, iupdwn, nupdwn, nudx, nelt
  USE electrons_module,         ONLY : distribute_c, collect_c
  USE efield_module,            ONLY : efield, epol, tefield, allocate_efield, &
                                       efield_update, ipolp, qmat, gqq, evalue,&
                                       berry_energy, pberryel, pberryion,      &
                                       efield2, epol2, tefield2,               &
                                       allocate_efield2, efield_update2,       &
                                       ipolp2, qmat2, gqq2, evalue2,           &
                                       berry_energy2, pberryel2, pberryion2
  USE ensemble_dft,             ONLY : tens, z0t, gibbsfe
  USE cg_module,                ONLY : tcg,  cg_update, c0old
  USE gvect,                    ONLY : ngm, ngm_g
  USE gvecs,                    ONLY : ngms
  USE smallbox_gvec,                    ONLY : ngb
  USE gvecw,                    ONLY : ngw, ngw_g
  USE gvect,       ONLY : gstart, mill, eigts1, eigts2, eigts3
  USE ions_base,                ONLY : na, nat, amass, nax, nsp, rcmax
  USE ions_base,                ONLY : ind_srt, ions_cofmass, ions_kinene, &
                                       ions_temp, ions_thermal_stress, &
                                       if_pos, extfor
  USE ions_base,                ONLY : ions_vrescal, fricp, greasp, &
                                       iforce, ndfrz, ions_shiftvar, ityp, &
                                       atm, ind_bck, cdm, cdms, ions_cofmsub
  USE cell_base,                ONLY : at, bg, ainv, frich, &
                                       greash, tpiba2, omega, alat, ibrav,  &
                                       celldm, h, hold, hnew, velh,         &
                                       wmass, press, iforceh, cell_force
  USE local_pseudo,             ONLY : allocate_local_pseudo
  USE io_global,                ONLY : stdout, ionode, ionode_id
  USE dener,                    ONLY : detot
  USE constants,                ONLY : pi, k_boltzmann_au, au_ps
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
                                       nhpbeg, nhpend,               &
                                       vnhp, xnhp0, xnhpm, xnhpp,    &
                                       atm2nhp, ions_nosevel, ions_noseupd,  &
                                       tempw, ions_nose_nrg, gkbt2nhp,       &
                                       ekin2nhp, anum2nhp
  USE electrons_nose,           ONLY : qne, ekincw, xnhe0, xnhep, xnhem,  &
                                       vnhe, electrons_nose_nrg,    &
                                       electrons_nose_shiftvar,           &
                                       electrons_nosevel, electrons_noseupd
  USE pres_ai_mod,              ONLY : P_ext, P_in, P_fin, pvar, volclu, &
                                       surfclu, Surf_t, abivol, abisur
  USE wavefunctions_module,     ONLY : c0_bgrp, cm_bgrp, phi_bgrp
  USE wannier_module,           ONLY : allocate_wannier
  USE cp_interfaces,            ONLY : printout_new, move_electrons, newinit
  USE cell_nose,                ONLY : xnhh0, xnhhm, xnhhp, vnhh, temph, &
                                       qnh, cell_nosevel, cell_noseupd,  &
                                       cell_nose_nrg, cell_nose_shiftvar
  USE cell_base,                ONLY : cell_kinene, cell_gamma, &
                                       cell_move, cell_hmove
  USE gvecw,                    ONLY : ecutwfc
  USE gvect,                    ONLY : ecutrho
  USE time_step,                ONLY : delt, tps, dt2,  twodelt
  USE cp_interfaces,            ONLY : cp_print_rho, nlfh, print_lambda, prefor, dotcsc
  USE cp_main_variables,        ONLY : acc, lambda, lambdam, lambdap, &
                                       ema0bg, sfac, eigr, iprint_stdout,  &
                                       irb, taub, eigrb, rhog, rhos, &
                                       rhor, bephi, becp_bgrp, nfi, descla, &
                                       drhor, drhog, bec_bgrp, dbec
  USE autopilot,                ONLY : event_step, event_index, &
                                       max_event_step, restart_p
  USE cell_base,                ONLY : s_to_r, r_to_s
  USE wannier_subroutines,      ONLY : wannier_startup, wf_closing_options, &
                                       ef_enthalpy
  USE cp_interfaces,            ONLY : writefile, eigs, strucf, phfacs
  USE cp_interfaces,            ONLY : ortho, elec_fakekine, calbec_bgrp, calbec, caldbec_bgrp
  USE constraints_module,       ONLY : check_constraint, remove_constr_force
  USE cp_autopilot,             ONLY : pilot
  USE ions_nose,                ONLY : ions_nose_allocate, ions_nose_shiftvar
  USE orthogonalize_base,       ONLY : updatc
  USE control_flags,            ONLY : force_pairing
  USE mp,                       ONLY : mp_bcast, mp_sum
  USE mp_global,                ONLY : root_bgrp, intra_bgrp_comm, np_ortho, &
                                       me_ortho, ortho_comm, &
                                       me_bgrp, inter_bgrp_comm, nbgrp, me_image
  USE ldaU_cp,                  ONLY : lda_plus_u, vupsi
  USE fft_base,                 ONLY : dfftp
  USE london_module,            ONLY : energy_london, force_london, stres_london
  USE input_parameters,         ONLY : tcpbo
  USE funct,                    ONLY : dft_is_hybrid, start_exx, exx_is_active
  !
  IMPLICIT NONE
  !
  ! ... input/output variables
  !
  REAL(DP), INTENT(OUT) :: tau_out(3,nat)
  REAL(DP), INTENT(OUT) :: fion_out(3,nat)
  REAL(DP), INTENT(OUT) :: etot_out
  !
  ! ... control variables
  !
  LOGICAL :: tfirst, tlast, tstop, tconv
  LOGICAL :: tprint, tfile, tstdout
    !  logical variable used to control printout
  !
  ! ... forces on ions
  !
  REAL(DP) :: maxfion, fion_tot(3)
  !
  ! ... work variables
  !
  REAL(DP) :: tempp, savee, saveh, savep, epot, epre, &
              enow, econs, econt, fccc, ccc, bigr, dt2bye
  REAL(DP) :: ekinc0, ekinp, ekinpr, ekinc
  REAL(DP) :: temps(nat)
  REAL(DP) :: ekinh, temphc, randy
  REAL(DP) :: delta_etot
  REAL(DP) :: ftmp, enb, enbi
  INTEGER  :: is, nacc, ia, j, iter, i, isa, ipos, iat, CYCLE_NOSE
  INTEGER  :: k, ii, l, m, iss
  REAL(DP) :: hgamma(3,3), temphh(3,3)
  REAL(DP) :: fcell(3,3)
  REAL(DP) :: deltaP, ekincf
  REAL(DP) :: stress_gpa(3,3), thstress(3,3), stress(3,3)
  !
  REAL(DP), ALLOCATABLE :: usrt_tau0(:,:), usrt_taup(:,:), usrt_fion(:,:)
    ! temporary array used to store unsorted positions and forces for
    ! constrained dynamics
  CHARACTER(LEN=3) :: labelw( nat )
    ! for force_pairing
  INTEGER   :: nspin_sub , i1, i2
    ! pmass contains masses in atomic Hartree units
  REAL(DP), ALLOCATABLE :: pmass(:)
  REAL(DP), ALLOCATABLE :: forceh(:,:)
  !
  CALL start_clock( 'cpr_total' )
  !
  etot_out = 0.D0
  enow     = 1.D9
  stress   = 0.0D0
  thstress   = 0.0D0
  !
  tfirst = .TRUE.
  tlast  = .FALSE.
  nacc   = 5
  !
  ALLOCATE ( pmass (nsp) )
  pmass(1:nsp) = amass(1:nsp) * amu_au
  nspin_sub = nspin
  IF( force_pairing ) nspin_sub = 1
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
  IF ( lda_plus_u ) ALLOCATE( forceh( 3, nat ) )
  !
  !
  !======================================================================
  !
  !           basic loop for molecular dynamics starts here
  !
  !======================================================================
  !
  main_loop: DO
     !
     CALL start_clock( 'main_loop' )
     CALL start_clock( 'cpr_md' )
     !
     dt2bye   = dt2 / emass
     nfi     = nfi + 1
     tlast   = ( nfi == nomore ) .OR. tlast
     tprint  = ( MOD( nfi, iprint ) == 0 ) .OR. tlast 
     tfile   = ( MOD( nfi, iprint ) == 0 )
     tstdout = ( MOD( nfi, iprint_stdout ) == 0 ) .OR. tlast
     !
     IF ( abivol ) THEN
        IF ( pvar ) THEN
           IF ( nfi .EQ. 1 ) THEN
              deltaP = (P_fin - P_in) / DBLE(nomore)
              P_ext = P_in
           ELSE
              P_ext = P_ext + deltaP
           END IF
        END IF
     END IF
     !
     IF ( ionode .AND. tstdout ) &
        WRITE( stdout, '(/," * Physical Quantities at step:",I6)' ) nfi
     !
     IF ( tnosee ) THEN
        fccc = 1.D0 / ( 1.D0 + 0.5D0 * delt * vnhe )
     ELSE IF ( tsde ) THEN
        fccc = 1.D0 
     ELSE
        fccc = 1.D0 / ( 1.D0 + frice )
     END IF
     !
     ! ... calculation of velocity of nose-hoover variables
     !
     IF ( tnosep ) THEN
        !
        CALL ions_nosevel( vnhp, xnhp0, xnhpm, delt, nhpcl, nhpdim )
        !
     END IF
     !
     IF ( tnosee ) THEN
        !
        CALL electrons_nosevel( vnhe, xnhe0, xnhem, delt )
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
     IF ( (okvan .or. nlcc_any ) .AND. (tfor .OR. thdyn .OR. tfirst) ) THEN
        !
        CALL initbox( tau0, alat, at, ainv, taub, irb )
        !
        CALL phbox( taub, iverbosity, eigrb )
        !
     END IF
     !
     IF ( tfor .OR. thdyn ) THEN
        !
        CALL phfacs( eigts1,eigts2,eigts3, eigr, mill, taus, dfftp%nr1,dfftp%nr2,dfftp%nr3, nat )
        !
        ! ... strucf calculates the structure factor sfac
        !
        CALL strucf( sfac, eigts1, eigts2, eigts3, mill, ngms )
        !
     END IF
     !
     IF ( thdyn ) THEN
        !
        CALL formf( tfirst, eself )
        !
     END IF
     !
     ! ... why this call ??? from Paolo Umari
     !
     IF ( tefield .or. tefield2 ) THEN
        !
        CALL calbec( 1, nsp, eigr, c0_bgrp, bec_bgrp ) ! ATTENZIONE  
        !
     END IF
     !
     ! Autopilot (Dynamic Rules) Implimentation    
     !
     call pilot(nfi)
     !
     IF ( ( tfor .OR. tfirst ) .AND. tefield ) CALL efield_update( eigr )
     IF ( ( tfor .OR. tfirst ) .AND. tefield2 ) CALL efield_update2( eigr )
     !
     ! ... pass ions information to plugins
     !
     CALL plugin_init_ions( tau0 )
     !
     IF ( lda_plus_u ) then
        ! forceh    ! Forces on ions due to Hubbard U 
        forceh=0.0d0
        ! vupsi     ! potentials on electrons due to Hubbard U
        vupsi=(0.0d0,0.0d0)
        CALL new_ns(c0_bgrp,eigr,vkb,vupsi,forceh)
        if ( mod(nfi,iprint).eq.0 ) call write_ns
     endif
     !
     !=======================================================================
     !
     !    electronic degrees of freedom are updated here
     !
     !=======================================================================
     !
     IF( force_pairing ) THEN
          c0_bgrp(:,iupdwn(2):nbsp)       =     c0_bgrp(:,1:nupdwn(2))
          cm_bgrp(:,iupdwn(2):nbsp)       =     cm_bgrp(:,1:nupdwn(2))
         phi_bgrp(:,iupdwn(2):nbsp)       =    phi_bgrp(:,1:nupdwn(2))
      lambda(:,:, 2) = lambda(:,:, 1)
     ENDIF
     !
     ! ... fake electronic kinetic energy
     !
     IF ( .NOT. tcg ) THEN
        !
        ekincf = 0.0d0

        CALL elec_fakekine( ekincf, ema0bg, emass, cm_bgrp, c0_bgrp, ngw, nbsp_bgrp, 1, delt )
        !
     END IF
     !
     CALL move_electrons( nfi, tfirst, tlast, bg(:,1), bg(:,2), bg(:,3), &
                          fion, c0_bgrp, cm_bgrp, phi_bgrp, &
                          enthal, enb, enbi, fccc, ccc, dt2bye, stress, .false. )
     !
     IF (lda_plus_u) fion = fion + forceh
     !
     ! DFT+D (Grimme) dispersion energy, forces (factor 0.5 converts to Ha/a.u.)
     !
     IF ( llondon ) THEN
        ALLOCATE( usrt_tau0( 3, nat ))
        usrt_tau0(:,:) = tau0(:,ind_bck(:))/alat
        delta_etot = 0.5_dp*energy_london (alat, nat,ityp,at,bg, usrt_tau0)
        etot = etot + delta_etot
        enthal=enthal+delta_etot
        IF ( tfor ) THEN
           ALLOCATE( usrt_fion( 3, nat ) )
           usrt_fion =  0.5_dp*force_london ( alat, nat,ityp, at,bg, usrt_tau0 )
           fion(:,:) = fion(:,:) + usrt_fion(:,ind_srt(:))
           DEALLOCATE (usrt_fion)
        END IF
        IF ( tpre ) stress = stress + 0.5_dp * stres_london ( alat , nat , &
                              ityp , at , bg , usrt_tau0 , omega )
        DEALLOCATE ( usrt_tau0 )
     END IF
     !
     IF ( tpre ) THEN
        !
        CALL nlfh( stress, bec_bgrp, dbec, lambda, descla )
        !
        CALL ions_thermal_stress( stress, thstress, pmass, omega, h, vels, nsp, na )
        !
        IF (tstdout) THEN
          WRITE(stdout,'(5X,"Pressure of Nuclei (GPa)",F20.5,I7)') &
              (thstress(1,1)+thstress(2,2)+thstress(3,3))/3.0_DP * au_gpa, nfi
          WRITE(stdout,'(5X,"Pressure Total (GPa)",F20.5,I7)') &
              (stress(1,1)+stress(2,2)+stress(3,3))/3.0_DP * au_gpa , nfi
        END IF
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
     ! BS : Initialization of additional cycles for the Nose thermostat ... 
     !
     IF (tnosep) CYCLE_NOSE=0
     !
444  IF ( tfor ) THEN
        !
        IF ( lwf ) CALL ef_force( fion, na, nsp, zv )
        !
        IF( textfor ) THEN 
           !
           FORALL( ia = 1:nat ) fion(:,ia) = fion(:,ia) + extfor(:,ia)
           !
           fion_tot(:) = SUM( fion(:,:), DIM = 2 ) / DBLE( nat )
           !
           FORALL( ia = 1:nat ) fion(:,ia) = fion(:,ia) - fion_tot(:)
           !
        END IF
        !
        IF ( remove_rigid_rot ) &
           CALL remove_tot_torque( nat, tau0, pmass(ityp(ind_srt(:))), fion )
        !
        IF ( lconstrain ) THEN
           !
           IF ( ionode ) THEN
              !
              ALLOCATE( usrt_tau0( 3, nat ) )
              ALLOCATE( usrt_taup( 3, nat ) )
              ALLOCATE( usrt_fion( 3, nat ) )
              !
              usrt_tau0(:,:) = tau0(:,ind_bck(:))
              usrt_fion(:,:) = fion(:,ind_bck(:))
              !
              ! ... we first remove the component of the force along the 
              ! ... constrain gradient (this constitutes the initial guess 
              ! ... for the lagrange multiplier)
              !
              CALL remove_constr_force( nat, usrt_tau0, if_pos, ityp, 1.D0, usrt_fion )
              !
              fion(:,:) = usrt_fion(:,ind_srt(:))
              !
           END IF
           !
           CALL mp_bcast( fion, ionode_id, intra_bgrp_comm )
           !
        END IF
        !
        !
        ! ... call void routine for user define/ plugin patches on external forces
        !
        CALL plugin_ext_forces()
        !
        !
        CALL ions_move( tausp, taus, tausm, iforce, pmass, fion, ainv, &
                        delt, na, nsp, fricp, hgamma, vels, tsdp, tnosep, &
                        fionm, vnhp, velsp, velsm, nhpcl, nhpdim, atm2nhp )
        !
        IF ( lconstrain ) THEN
           !
           ! ... constraints are imposed here
           !
           IF ( ionode ) THEN
              !
              CALL s_to_r( tausp, taup, na, nsp, hnew )
              !
              usrt_taup(:,:) = taup(:,ind_bck(:))
              !
              CALL check_constraint( nat, usrt_taup, usrt_tau0, usrt_fion, &
                                     if_pos, ityp, 1.D0, delt, amu_au )
              !
              taup(:,:) = usrt_taup(:,ind_srt(:))
              fion(:,:) = usrt_fion(:,ind_srt(:))
              !
              DEALLOCATE( usrt_tau0, usrt_taup, usrt_fion )
              !
           END IF
           !
           CALL mp_bcast( taup, ionode_id, intra_bgrp_comm )
           CALL mp_bcast( fion, ionode_id, intra_bgrp_comm )
           !
           CALL r_to_s( taup, tausp, na, nsp, ainv )
           !
        END IF
        !
        CALL ions_cofmass( tausp, pmass, na, nsp, cdm )
        !
        IF ( ndfrz == 0 ) &
           CALL ions_cofmsub( tausp, iforce, nat, cdm, cdms )
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
     !
     IF ( tfor .OR. thdyn ) THEN
        !
        IF ( .NOT.tnosep .OR. CYCLE_NOSE.EQ.0 ) THEN
          !
          IF ( thdyn ) THEN
             !
             hold = h
             h    = hnew
             !
             IF( nbgrp > 1 ) THEN
                CALL mp_bcast( h, 0, inter_bgrp_comm )
             END IF
             !
             CALL newinit( h, iverbosity )
             !
             CALL newnlinit()
             !
          ELSE
             !
             hold = h
             !
          END IF
          !
        END IF
        !
        ! ... phfac calculates eigr
        !
        CALL phfacs( eigts1,eigts2,eigts3, eigr, mill, tausp, dfftp%nr1,dfftp%nr2,dfftp%nr3, nat ) 
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
     ! In case of tnosep = .true., the orthonormality is done only with the most updated 
     ! atomic coordinates coming out of the CYCLE_NOSE loop
     !
     IF ( .NOT.tnosep .OR. CYCLE_NOSE.EQ.2 ) THEN
       !
       IF ( .NOT. tcg ) THEN
         !
         IF ( tortho ) THEN
           !
           CALL ortho( eigr, cm_bgrp, phi_bgrp, lambda, descla, bigr, iter, ccc, bephi, becp_bgrp )
           !
         ELSE
           !
           CALL gram_bgrp( vkb, bec_bgrp, nkb, cm_bgrp, ngw )
           !
           IF ( iverbosity > 2 ) CALL dotcsc( eigr, cm_bgrp, ngw, nbsp_bgrp )
           !
         END IF
         !
         !  correction to displacement of ions
         !
         IF ( iverbosity > 1 ) CALL print_lambda( lambda, descla, nbsp, 9, 1.D0 )
         !
         IF ( tortho ) THEN
           CALL updatc( ccc, lambda, phi_bgrp, bephi, becp_bgrp, bec_bgrp, cm_bgrp, descla )
         END IF
         !
         IF( force_pairing ) THEN
           c0_bgrp(:,iupdwn(2):nbsp)       =     c0_bgrp(:,1:nupdwn(2))
           cm_bgrp(:,iupdwn(2):nbsp)       =     cm_bgrp(:,1:nupdwn(2))
           phi_bgrp(:,iupdwn(2):nbsp)       =    phi_bgrp(:,1:nupdwn(2))
           lambda(:,:, 2) = lambda(:,:, 1)
         ENDIF
         !
         CALL calbec_bgrp( nvb+1, nsp, eigr, cm_bgrp, bec_bgrp )
         !
         IF ( tpre ) THEN
           CALL caldbec_bgrp( eigr, cm_bgrp, dbec, descla )
         END IF
         !
         IF ( iverbosity > 1 ) CALL dotcsc( eigr, cm_bgrp, ngw, nbsp_bgrp )
         !
       END IF
       !
     END IF !(.NOT.tnosep.OR.CYCLE_NOSE.EQ.2)
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
     ! ... ionic kinetic energy and temperature
     !
     IF ( tfor ) THEN
        !
        CALL ions_vel( vels, tausp, tausm, na, nsp, delt )
        !
        CALL ions_kinene( ekinp, vels, na, nsp, hold, pmass )
        !
        CALL ions_temp( tempp, temps, ekinpr, vels, na, nsp, &
                        hold, pmass, ndega, nhpdim, atm2nhp, ekin2nhp )
        !
     END IF
     !
     ! ... fake electronic kinetic energy
     !
     IF ( .NOT. tcg ) THEN
        !
        CALL elec_fakekine( ekinc0, ema0bg, emass, c0_bgrp, cm_bgrp, ngw, nbsp_bgrp, 1, delt )
        !
        ekinc0 = (ekinc0 + ekincf)*0.5d0
        !
        ekinc = ekinc0
        !
     END IF
     !
     ! ... fake cell-parameters kinetic energy
     !
     ekinh = 0.D0
     !
     IF ( thdyn ) THEN
        !
        CALL cell_kinene( ekinh, temphh, velh )
        !
     END IF
     !
     IF ( COUNT( iforceh == 1 ) > 0 ) THEN
        !
        temphc = 2.D0 / k_boltzmann_au * ekinh / DBLE( COUNT( iforceh == 1 ) )
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
                                      nhpcl, nhpdim, nhpbeg, nhpend )
     !
     IF ( tnosee ) CALL electrons_noseupd( xnhep, xnhe0, xnhem, &
                                           delt, qne, ekinc, ekincw, vnhe )
     !
     IF ( tnoseh ) CALL cell_noseupd( xnhhp, xnhh0, xnhhm, &
                                      delt, qnh, temphh, temph, vnhh )
     !
     !=================================================================
     ! BS : Additional cycles for the Nose thermostat ... 
     IF(tnosep) CYCLE_NOSE=CYCLE_NOSE+1
     IF(tnosep .AND. (CYCLE_NOSE .LE. 2) ) GO TO 444 
     !=================================================================
     ! 
     ! ... warning:  thdyn and tcp/tcap are not compatible yet!!!
     !
     IF ( tcp .AND. tfor .AND. .NOT.thdyn ) THEN
        !
        IF ( tempp > (tempw+tolp) .OR. &
             tempp < (tempw-tolp) .AND. tempp /= 0.D0 ) THEN
           !
           CALL  ions_vrescal( tcap, tempw, tempp, taup, &
                               tau0, taum, na, nsp, fion, iforce, pmass, delt )
           CALL r_to_s( taup, tausp, na, nsp, ainv ) 
           !
        END IF
        !
     END IF
     !
     IF ( tprint ) THEN
        !
        IF( tortho ) THEN
           !
           IF( force_pairing )  THEN
              lambda(:, :, 2) =  lambda(:, :, 1)
              lambdap(:, :, 2) = lambdap(:, :, 1)
              WRITE( stdout, '("Occupations in CPR:")' )
              WRITE( stdout, '(10F9.6)' ) ( f(i), i = 1, nbspx )  
           END IF
           !
           CALL eigs( nfi, lambdap, lambda, descla )
           !
        ELSE
           !
           WRITE( stdout, '("NOTE: eigenvalues are not computed without ortho")' )
           !
        END IF
        !
     END IF
     !
     IF ( lwf ) CALL ef_enthalpy( enthal, tau0 )
     !
     IF ( tens .AND. tprint ) THEN
        !
        WRITE( stdout, '("Occupations  :")' )
        WRITE( stdout, '(10F9.6)' ) ( f(i), i = 1, nbsp )
        !
     END IF
     !
     epot = eht + epseu + exc
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
     if (abivol) etot = etot - P_ext*volclu
     if (abisur) etot = etot - Surf_t*surfclu
     !
     IF ( tstdout) CALL spinsq ( c0_bgrp, bec_bgrp, rhor )
     !
     CALL printout_new( nfi, tfirst, tfile, tprint, tps, hold, stress, &
                        tau0, vels, fion, ekinc, temphc, tempp, temps, etot, &
                        enthal, econs, econt, vnhh, xnhh0, vnhp, xnhp0, vnhe, xnhe0, atot, &
                        ekin, epot, tprnfor, tpre, tstdout )
     !
     if (abivol) etot = etot + P_ext*volclu
     if (abisur) etot = etot + Surf_t*surfclu
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
        IF( nbgrp > 1 ) THEN
           CALL mp_bcast( tau0, 0, inter_bgrp_comm )
           CALL mp_bcast( taus, 0, inter_bgrp_comm )
           CALL mp_bcast( vels, 0, inter_bgrp_comm )
           CALL mp_bcast( xnhp0, 0, inter_bgrp_comm )
           CALL mp_bcast( xnhe0, 0, inter_bgrp_comm )
           CALL mp_bcast( xnhh0, 0, inter_bgrp_comm )
        END IF
        !
     END IF
     !
     ekincm = ekinc0
     !  
     ! ... cm=c(t+dt) c0=c(t)
     !
     IF( .NOT. tcg ) THEN
        !
        CALL dswap( 2*SIZE( c0_bgrp ), c0_bgrp, 1, cm_bgrp, 1 )
        !
     ELSE
        !
        CALL cg_update( tfirst, nfi, c0_bgrp )
        !
        IF ( tfor .AND. .NOT. tens .AND. tprint ) THEN
           !
           ! ... in this case optimize c0 and lambda for smooth
           ! ... restart with CP
           !
           IF ( okvan .or. nlcc_any ) THEN
              CALL initbox( tau0, alat, at, ainv, taub, irb )
              CALL phbox( taub, iverbosity, eigrb ) 
           END IF
           CALL r_to_s( tau0, taus, na, nsp, ainv )
           CALL phfacs( eigts1,eigts2,eigts3, eigr, mill, taus, dfftp%nr1,dfftp%nr2,dfftp%nr3, nat )
           CALL strucf( sfac, eigts1, eigts2, eigts3, mill, ngms )
           !
           IF ( thdyn )    CALL formf( tfirst, eself )
           IF ( tefield )  CALL efield_update( eigr )
           IF ( tefield2 ) CALL efield_update2( eigr )
           !
           CALL plugin_init_ions( tau0 )
           !
           lambdam = lambda
           !
           CALL move_electrons( nfi, tfirst, tlast, bg(:,1), bg(:,2), bg(:,3),&
                                fion, c0_bgrp, cm_bgrp, phi_bgrp, enthal, enb,&
                                enbi, fccc, ccc, dt2bye, stress,.true. )
           !
        END IF
        !
     END IF
     !
     ! ... now:  cm=c(t) c0=c(t+dt)
     !
     tfirst = .FALSE.
     !
     CALL stop_clock( 'main_loop' )
     !
     ! ... write on file ndw each isave
     !
     IF ( ( MOD( nfi, isave ) == 0 ) .AND. ( nfi < nomore ) ) THEN
        !
        IF ( tcg ) THEN
          !
          CALL writefile( h, hold ,nfi, c0_bgrp, c0old, taus, tausm,  &
                          vels, velsm, acc, lambda, lambdam, descla, xnhe0, xnhem,     &
                          vnhe, xnhp0, xnhpm, vnhp, nhpcl,nhpdim,ekincm, xnhh0,&
                          xnhhm, vnhh, velh, fion, tps, z0t, f, rhor )
           !
        ELSE
           !
           CALL writefile( h, hold, nfi, c0_bgrp, cm_bgrp, taus,  &
                           tausm, vels, velsm, acc,  lambda, lambdam, descla, xnhe0,   &
                           xnhem, vnhe, xnhp0, xnhpm, vnhp, nhpcl, nhpdim, ekincm,&
                           xnhh0, xnhhm, vnhh, velh, fion, tps, z0t, f, rhor )
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
     delta_etot = ABS( epre - enow )
     !
     !exx_wf related
     !The following criteria is used to turn on exact exchange calculation when
     !GGA energy is converged up to 100 times of the input etot convergence thereshold  
     !
     IF( .NOT.exx_is_active().AND.dft_is_hybrid().AND.tconvthrs%active ) THEN
       !
       IF(delta_etot.LT.tconvthrs%derho*1.E+2_DP) THEN
         !
         WRITE(stdout,'(/,3X,"Exact Exchange is turned on ...")')
         ! 
         CALL start_exx()
         !
       END IF
       !
     END IF
     !
     tstop = check_stop_now() .OR. tlast
     !
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
        tlast = .TRUE.
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
     IF ( lwf ) &
        CALL wf_closing_options( nfi, c0_bgrp, cm_bgrp, bec_bgrp, eigr, eigrb,&
                                 taub, irb, ibrav, bg(:,1), bg(:,2), bg(:,3), &
                                 taus, tausm, vels, &
                                 velsm, acc, lambda, lambdam, descla, xnhe0, xnhem,  &
                                 vnhe, xnhp0, xnhpm, vnhp, nhpcl, nhpdim,    &
                                 ekincm, xnhh0, xnhhm, vnhh, velh, ecutrho,  &
                                 ecutwfc,delt,celldm, fion, tps, z0t, f, rhor )
     !
     IF ( tstop ) EXIT main_loop
     !
     CALL stop_clock( 'cpr_md' )
     !
  END DO main_loop
  !
  !===================== end of main loop of molecular dynamics ===============
  !
  DEALLOCATE ( pmass )
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
        ipos = ind_srt( isa )
        tau_out(:,ipos) = tau0(:,isa)
        fion_out(:,ipos) = fion(:,isa)
        !
     END DO
     !
  END DO
  !
  conv_elec = .TRUE.
  !
  IF ( tcg ) cm_bgrp = c0old
  !
  CALL writefile( h, hold, nfi, c0_bgrp, cm_bgrp, taus, tausm, &
                  vels, velsm, acc, lambda, lambdam, descla, xnhe0, xnhem, vnhe,    &
                  xnhp0, xnhpm, vnhp, nhpcl,nhpdim,ekincm, xnhh0, xnhhm,    &
                  vnhh, velh, fion, tps, z0t, f, rhor )
  !
  IF( iverbosity > 1 ) CALL print_lambda( lambda, descla, nbsp, nbsp, 1.D0 )
  !
  IF (lda_plus_u) DEALLOCATE( forceh )

  !
  CALL stop_clock( 'cpr_total' ) ! BS
  !
  RETURN
  !
END SUBROUTINE cprmain
!
!----------------------------------------------------------------------------
SUBROUTINE terminate_run()
  !----------------------------------------------------------------------------
  !
  USE io_global,         ONLY : stdout, ionode
  USE control_flags,     ONLY : ts_vdw, thdyn, tortho
  USE cg_module,         ONLY : tcg, print_clock_tcg
  USE ldaU_cp,           ONLY : lda_plus_u
  USE mp,                ONLY : mp_report
  USE control_flags,     ONLY : lwf, lwfpbe0nscf
  USE tsvdw_module,      ONLY : tsvdw_finalize
  USE exx_module,        ONLY : exx_finalize
  USE funct,             ONLY : dft_is_hybrid, exx_is_active
  !
  IMPLICIT NONE
  !
  ! ...  print statistics
  !
  CALL printacc()
  !
  !==============================================================
  WRITE( stdout, '(/5x,"Called by MAIN_LOOP:")' )
  CALL print_clock( 'total_time' )
  CALL print_clock( 'initialize' )
  CALL print_clock( 'main_loop' )
  CALL print_clock( 'cpr_total' )
  !==============================================================
  WRITE( stdout, '(/5x,"Called by INIT_RUN:")' )
  CALL print_clock( 'init_readfile' )
  IF( lwf ) CALL print_clock( 'wf_start' )
  !==============================================================
  WRITE( stdout, '(/5x,"Called by CPR:")' )
  CALL print_clock( 'cpr_md' )
  CALL print_clock( 'move_electrons' )
  CALL print_clock( 'move_ion' )
  IF( lwf ) CALL print_clock( 'wf_close_opt' ) ! wf_close_opt = wf_1 + wf_2
  !==============================================================
  IF( lwf ) THEN
    WRITE( stdout, '(/5x,"Called by WANNIER_MODULES:")' )
    CALL print_clock( 'wf_start' )
    CALL print_clock( 'wf_init' )
    CALL print_clock( 'wf_close_opt' ) ! wf_close_opt = wf_1 + wf_2
    CALL print_clock( 'wf_1' )
    CALL print_clock( 'wf_2' )
    CALL print_clock( 'ddyn_u' )
    CALL print_clock( 'ortho_u' )
  END IF
  !==============================================================
  !exx_wf related
  IF ( dft_is_hybrid().AND.exx_is_active() ) THEN
    !
    WRITE( stdout, '(/5x,"Called by EXACT_EXCHANGE:")' )
    CALL print_clock('exact_exchange')   ! total time for exx
    CALL print_clock('self_v')           ! calculation self potential
    CALL print_clock('getpairv')         ! calculation pair potential
    CALL print_clock('exx_gs_setup')     ! calculation
    CALL print_clock('exx_pairs')        ! calculation
    CALL print_clock('r_orbital')        ! communication 
    CALL print_clock('totalenergy')      ! communication 
    CALL print_clock('vl2vg')            ! communication
    CALL print_clock('send_psi')         ! communication
    CALL print_clock('sendv')            ! communication
    CALL print_clock('send_psi_barrier') ! communication
    CALL print_clock('send_psi_wait')    ! communication
    !CALL print_clock('sendv_barrier')    ! communication
    CALL print_clock('getvofr')
    CALL print_clock('getvofr_qlm')
    CALL print_clock('getvofr_bound')
    CALL print_clock('getvofr_geterho')
    CALL print_clock('getvofr_hpotcg')
    !CALL print_clock('hpotcg')
    !CALL print_clock('getqlm')
    !CALL print_clock('exx_bound')
    !CALL print_clock('lapmv')
    CALL print_clock('exx_cell_derv')
    !
    CALL exx_finalize() ! deallocate all arrays
    !
  END IF
  !==============================================================
  IF (thdyn) THEN 
    WRITE( stdout, '(/5x,"Called by CELL_DYNAMICS:")' )
    CALL print_clock( 'formf' )
  END IF

  WRITE( stdout, '(/5x,"Called by move_electrons:")' )
  CALL print_clock( 'rhoofr' )
  CALL print_clock( 'vofrho' )
  CALL print_clock( 'dforce' )
  CALL print_clock( 'calphi' )
  CALL print_clock( 'newd' )
  CALL print_clock( 'nlfl' )

  IF (lda_plus_u) WRITE( stdout, '(/5x,"Called by new_ns:")' )
  CALL print_clock( 'new_ns:forc' )
  CALL print_clock( 'projwfc_hub' )
  CALL print_clock( 'dndtau' )

  IF (tortho) WRITE( stdout, '(/5x,"Called by ortho:")' )
  IF (tortho) THEN
    CALL print_clock( 'ortho_iter' )
    CALL print_clock( 'rsg' )
    CALL print_clock( 'rhoset' )
    CALL print_clock( 'sigset' )
    CALL print_clock( 'tauset' )
    CALL print_clock( 'ortho' )
    CALL print_clock( 'updatc' )
  ELSE
    CALL print_clock( 'gram' )
  END IF

  WRITE( stdout, '(/5x,"Small boxes:")' )
  CALL print_clock( 'rhov' )
  CALL print_clock( 'fftb' )
  CALL print_clock( 'set_cc' )
  CALL print_clock( 'forcecc' )

  WRITE( stdout, '(/5x,"Low-level routines:")' )
  CALL print_clock( 'prefor' )
  CALL print_clock( 'nlfq' )
  CALL print_clock( 'nlsm1' )
  CALL print_clock( 'nlsm2' )
  CALL print_clock( 'fft' )
  CALL print_clock( 'ffts' )
  CALL print_clock( 'fftw' )
  CALL print_clock( 'fft_scatter' )
  CALL print_clock( 'betagx' )
  CALL print_clock( 'qradx' )
  CALL print_clock( 'tmp_clk1' )
  CALL print_clock( 'tmp_clk2' )
  CALL print_clock( 'tmp_clk3' )
  CALL print_clock( 'gram' )
  CALL print_clock( 'nlinit' )
  CALL print_clock( 'init_dim' )
  CALL print_clock( 'newnlinit' )
  CALL print_clock( 'from_scratch' )
  CALL print_clock( 'from_restart' )
  CALL print_clock( 'new_ns' )
  CALL print_clock( 'strucf' )
  CALL print_clock( 'calbec' )
!==============================================================
  IF (ts_vdw) THEN
    WRITE( stdout, '(/5x,"Called by tsvdw:")' )
    CALL print_clock( 'ts_vdw' )
    CALL print_clock( 'tsvdw_pair' )
    CALL print_clock( 'tsvdw_rhotot' )
    CALL print_clock( 'tsvdw_screen' )
    CALL print_clock( 'tsvdw_veff' )
    CALL print_clock( 'tsvdw_dveff' )
    CALL print_clock( 'tsvdw_energy' )
    CALL print_clock( 'tsvdw_wfforce' )
    !
    CALL tsvdw_finalize()
  END IF
  !
  IF (tcg) call print_clock_tcg()
  !
  CALL print_clock( 'ALLTOALL' )
  !
  CALL plugin_clock()
  !
  CALL mp_report()
  !
END SUBROUTINE terminate_run
