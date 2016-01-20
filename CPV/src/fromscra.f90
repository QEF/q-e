!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE from_scratch( )
    !
    USE kinds,                ONLY : DP
    USE control_flags,        ONLY : tranp, trane, iverbosity, tpre, tv0rd, &
                                     tfor, thdyn, &
                                     lwf, tprnfor, tortho, amprp, ampre,  &
                                     tsde, ortho_eps, ortho_max, &
                                     force_pairing
    USE ions_positions,       ONLY : taus, tau0, tausm, vels, velsm, fion, fionm
    USE ions_base,            ONLY : na, nsp, randpos, zv, ions_vel, vel_srt
    USE ions_base,            ONLY : cdmi, nat, iforce
    USE ions_nose,            ONLY : xnhp0, xnhpm, vnhp
    USE cell_base,            ONLY : ainv, h, s_to_r, ibrav, omega, press, &
                                     hold, r_to_s, deth, wmass, iforceh,   &
                                     cell_force, velh, at, alat
    USE cell_nose,            ONLY : xnhh0, xnhhm, vnhh
    USE electrons_nose,       ONLY : xnhe0, xnhem, vnhe
    use electrons_base,       ONLY : nbsp, f, nspin, nupdwn, iupdwn, nbsp_bgrp, nbspx_bgrp, nbspx
    USE electrons_module,     ONLY : occn_info, distribute_c, collect_c, distribute_b, collect_b
    USE energies,             ONLY : entropy, eself, enl, ekin, enthal, etot, ekincm
    USE energies,             ONLY : dft_energy_type, debug_energies
    USE dener,                ONLY : denl, denl6, dekin6, detot
    USE uspp,                 ONLY : vkb, becsum, deeq, nkb, okvan, nlcc_any
    USE io_global,            ONLY : stdout, ionode
    USE core,                 ONLY : rhoc
    USE gvecw,                ONLY : ngw
    USE gvecs,                ONLY : ngms
    USE gvect,                ONLY : ngm, gg
    USE gvect,                ONLY : gstart, mill, eigts1, eigts2, eigts3
    USE uspp_param,           ONLY : nvb
    USE cp_electronic_mass,   ONLY : emass
    USE efield_module,        ONLY : tefield, efield_berry_setup, berry_energy, &
                                     tefield2, efield_berry_setup2, berry_energy2
    USE cg_module,            ONLY : tcg
    USE ensemble_dft,         ONLY : tens, compute_entropy
    USE cp_interfaces,        ONLY : runcp_uspp, runcp_uspp_force_pairing, &
                                     strucf, phfacs, nlfh, vofrho, nlfl_bgrp, prefor
    USE cp_interfaces,        ONLY : rhoofr, ortho, wave_rand_init, elec_fakekine
    USE cp_interfaces,        ONLY : compute_stress, dotcsc, calbec_bgrp, caldbec_bgrp
    USE cp_interfaces,        ONLY : print_lambda, nlfq_bgrp, setval_lambda
    USE printout_base,        ONLY : printout_pos
    USE orthogonalize_base,   ONLY : updatc, calphi_bgrp
    USE wave_base,            ONLY : wave_steepest
    USE wavefunctions_module, ONLY : c0_bgrp, cm_bgrp, phi_bgrp
    USE fft_base,             ONLY : dfftp
    USE time_step,            ONLY : delt
    USE cp_main_variables,    ONLY : descla, bephi, becp_bgrp, nfi, &
                                     sfac, eigr, taub, irb, eigrb, bec_bgrp, &
                                     lambda, lambdam, lambdap, ema0bg, rhog, rhor, rhos, &
                                     vpot, ht0, edft, becdr_bgrp, dbec, drhor, drhog
    USE mp_global,            ONLY : np_ortho, me_ortho, ortho_comm, inter_bgrp_comm, nbgrp
    USE mp,                   ONLY : mp_sum
    USE matrix_inversion
    !
    IMPLICIT NONE
    !
    REAL(DP),    ALLOCATABLE :: emadt2(:), emaver(:)
    REAL(DP)                 :: verl1, verl2
    REAL(DP)                 :: bigr, dum
    INTEGER                  :: i, j, iter, iss, ierr, nspin_wfc
    LOGICAL                  :: tlast = .FALSE.
    REAL(DP)                 :: gam(1,1,1)
    REAL(DP)                 :: fcell(3,3), ccc, enb, enbi, fccc
    LOGICAL                  :: ttforce
    LOGICAL                  :: tstress
    LOGICAL, PARAMETER       :: ttprint = .TRUE.
    REAL(DP)                 :: ei_unp  
    REAL(DP)                 :: dt2bye
    INTEGER                  :: n_spin_start 
    LOGICAL                  :: tfirst = .TRUE.
    REAL(DP)                 :: stress(3,3)
    INTEGER                  :: i1, i2 
    !
    ! ... Subroutine body
    !
    CALL start_clock( 'from_scratch' )
    !
    nfi = 0
    !
    ttforce = tfor  .or. tprnfor
    tstress = thdyn .or. tpre
    !
    stress = 0.0d0
    !
    IF( tsde ) THEN
       fccc = 1.0d0
    ELSE
       fccc = 0.5d0
    END IF
    !
    dt2bye = delt * delt / emass
    !
    IF( ANY( tranp( 1:nsp ) ) ) THEN
       !
       CALL invmat( 3, h, ainv, deth )
       !
       CALL randpos( taus, na, nsp, tranp, amprp, ainv, iforce )
       !
       CALL s_to_r( taus, tau0, na, nsp, h )
       !
    END IF
    !
    CALL phfacs( eigts1, eigts2, eigts3, eigr, mill, taus, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat )
    !
    CALL strucf( sfac, eigts1, eigts2, eigts3, mill, ngms )
    !     
    IF ( okvan .OR. nlcc_any ) THEN
       CALL initbox ( tau0, alat, at, ainv, taub, irb )
       CALL phbox( taub, iverbosity, eigrb )
    END IF
    !
    !     pass ions informations to plugins
    !
    CALL plugin_init_ions( tau0 )
    !
    !     wfc initialization with random numbers
    !     
    CALL wave_rand_init( cm_bgrp )
    !
    IF ( ionode ) &
       WRITE( stdout, fmt = '(//,3X, "Wave Initialization: random initial wave-functions" )' )
    !
    ! ... prefor calculates vkb (used by gram)
    !
    CALL prefor( eigr, vkb )
    !
    nspin_wfc = nspin
    IF( force_pairing ) nspin_wfc = 1

    CALL gram_bgrp( vkb, bec_bgrp, nkb, cm_bgrp, ngw )

    IF( force_pairing ) cm_bgrp(:,iupdwn(2):iupdwn(2)+nupdwn(2)-1) = cm_bgrp(:,1:nupdwn(2))
    !
    if( iverbosity > 1 ) CALL dotcsc( eigr, cm_bgrp, ngw, nbsp )
    !
    ! ... initialize bands
    !
    CALL occn_info( f )
    !
    hold = h
    velh = 0.0d0
    fion = 0.0d0
    !
    IF ( tv0rd .AND. tfor ) THEN
       !
       ! ... vel_srt=starting velocities, read from input, are brough to
       ! ... scaled axis and copied into array vels. Since velocites are
       ! ... not actually used by the Verlet algorithm, we set tau(t-dt)
       ! ... to tausm=tau(t)-v*delta t so that the Verlet algorithm will 
       ! ... start with the correct velocity
       !
       CALL r_to_s( vel_srt, vels, na, nsp, ainv )
       tausm(:,:) =  taus(:,:) - vels(:,:)*delt
       velsm(:,:) =  vels(:,:)
    ELSE
       vels = 0.D0
       tausm = taus
    END IF
    !
    ! ... compute local form factors
    !
    CALL formf( tfirst, eself )
    !
    edft%eself = eself
    IF( tefield ) THEN
      CALL efield_berry_setup( eigr, tau0 )
    END IF
    IF( tefield2 ) THEN
      CALL efield_berry_setup2( eigr, tau0 )
    END IF
    !
    IF( .NOT. tcg ) THEN
       !
       CALL calbec_bgrp ( 1, nsp, eigr, cm_bgrp, bec_bgrp )
       !
       if ( tstress ) CALL caldbec_bgrp( eigr, cm_bgrp, dbec, descla )
       !
       CALL rhoofr( nfi, cm_bgrp, irb, eigrb, bec_bgrp, dbec, becsum, rhor, drhor, rhog, drhog, rhos, enl, denl, ekin, dekin6 )
       !
       edft%enl  = enl
       edft%ekin = ekin
       !
    END IF
    !
    !     put core charge (if present) in rhoc(r)
    !
    if ( nlcc_any ) CALL set_cc( irb, eigrb, rhoc )
    !
    IF( .NOT. tcg ) THEN
   
      IF( tens ) THEN
        CALL compute_entropy( entropy, f(1), nspin )
        entropy = entropy * nbsp
      END IF
      !
      vpot = rhor
      !
      CALL vofrho( nfi, vpot, drhor, rhog, drhog, rhos, rhoc, tfirst, tlast, &
     &  eigts1, eigts2, eigts3, irb, eigrb, sfac, tau0, fion )

      IF( tefield ) THEN
        CALL berry_energy( enb, enbi, bec_bgrp, cm_bgrp, fion ) 
        etot = etot + enb + enbi
      END IF
      IF( tefield2 ) THEN
        CALL berry_energy2( enb, enbi, bec_bgrp, cm_bgrp, fion )
        etot = etot + enb + enbi
      END IF

      CALL compute_stress( stress, detot, h, omega )

      if( iverbosity > 1 ) &
             CALL printout_pos( stdout, fion, nat, head = ' fion ' )

      CALL newd( vpot, irb, eigrb, becsum, fion )
      !
      IF( force_pairing ) THEN
         !
         CALL runcp_uspp_force_pairing( nfi, fccc, ccc, ema0bg, dt2bye, rhos,&
                    bec_bgrp, cm_bgrp, c0_bgrp, ei_unp, fromscra = .TRUE. )
         !
         CALL setval_lambda( lambda(:,:,2), nupdwn(1), nupdwn(1), 0.d0, descla(1) )
         !
      ELSE
         !
         CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec_bgrp, cm_bgrp, c0_bgrp, fromscra = .TRUE. )
         !
      ENDIF
      !
      !     nlfq needs deeq bec
      !
      IF( ttforce ) THEN
         CALL nlfq_bgrp( cm_bgrp, eigr, bec_bgrp, becdr_bgrp, fion )
      END IF
      !
      !     calphi calculates phi
      !     the electron mass rises with g**2
      !
      CALL calphi_bgrp( cm_bgrp, ngw, bec_bgrp, nkb, vkb, phi_bgrp, nbspx_bgrp, ema0bg )
      !
      IF( force_pairing ) &
         &   phi_bgrp( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) =    phi_bgrp( :, 1:nupdwn(2))

      if( tortho ) then
         CALL ortho( eigr, c0_bgrp, phi_bgrp, lambda, descla, bigr, iter, ccc, bephi, becp_bgrp )
      else
         CALL gram_bgrp( vkb, bec_bgrp, nkb, c0_bgrp, ngw )
      endif
      !
      IF ( ttforce ) THEN
         CALL nlfl_bgrp( bec_bgrp, becdr_bgrp, lambda, descla, fion )
      END IF

      if ( iverbosity > 1 ) CALL print_lambda( lambda, descla, nbsp, 9, ccc )

      !
      if ( tstress ) CALL nlfh( stress, bec_bgrp, dbec, lambda, descla )
      !
      IF ( tortho ) THEN
         CALL updatc( ccc, lambda, phi_bgrp, bephi, becp_bgrp, bec_bgrp, c0_bgrp, descla )
      END IF
      !
      IF( force_pairing ) THEN
         !
         c0_bgrp ( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) = c0_bgrp( :, 1:nupdwn(2))
         phi_bgrp( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) = phi_bgrp( :, 1:nupdwn(2))
         lambda(:,:,2) = lambda(:,:,1)
         !
      ENDIF
      !
      !
      CALL calbec_bgrp ( nvb+1, nsp, eigr, c0_bgrp, bec_bgrp )
      !
      if ( tstress ) CALL caldbec_bgrp( eigr, cm_bgrp, dbec, descla )

      if ( iverbosity > 1 ) CALL dotcsc( eigr, c0_bgrp, ngw, nbsp_bgrp )
      !
      xnhp0 = 0.0d0
      xnhpm = 0.0d0
      vnhp  = 0.0d0
      fionm = 0.0d0
      !
      CALL ions_vel( vels, taus, tausm, na, nsp, delt )
      !
      xnhh0(:,:) = 0.0d0
      xnhhm(:,:) = 0.0d0
      vnhh (:,:) = 0.0d0
      velh (:,:) = ( h(:,:) - hold(:,:) ) / delt
      !
      CALL elec_fakekine( ekincm, ema0bg, emass, c0_bgrp, cm_bgrp, ngw, nbsp_bgrp, 1, delt )

      xnhe0 = 0.0d0
      xnhem = 0.0d0
      vnhe  = 0.0d0

      lambdam = lambda
      !
    ELSE 
       !
       c0_bgrp = cm_bgrp
       !
    END IF
    !
    CALL stop_clock( 'from_scratch' )
    !
    RETURN
    !
END SUBROUTINE from_scratch
