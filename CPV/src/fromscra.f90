!
! Copyright (C) 2002-2024 Quantum ESPRESSO Foundation
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
                                     tsde, force_pairing, tcap
    USE input_parameters,     ONLY : startingwfc
    USE ions_positions,       ONLY : taus, tau0, tausm, vels, velsm, fion, fionm, &
                                     taum 
    USE ions_base,            ONLY : na, nsp, randpos, zv, ions_vel, vel, ityp, &
                                     amass, randvel
    USE ions_base,            ONLY : cdmi, nat, iforce
    USE ions_nose,            ONLY : xnhp0, xnhpm, vnhp, tempw
    USE cell_base,            ONLY : ainv, h, s_to_r, ibrav, omega, press, &
                                     hold, r_to_s, deth, wmass, iforceh,   &
                                     cell_force, velh, at, alat
    USE cell_nose,            ONLY : xnhh0, xnhhm, vnhh
    USE electrons_nose,       ONLY : xnhe0, xnhem, vnhe
    use electrons_base,       ONLY : nbsp, f, nspin, nupdwn, iupdwn, nbsp_bgrp, nbspx_bgrp, nbspx, nudx
    USE electrons_module,     ONLY : occn_info, distribute_c, collect_c, distribute_b, collect_b
    USE energies,             ONLY : entropy, eself, enl, ekin, enthal, etot, ekincm
    USE energies,             ONLY : dft_energy_type, debug_energies
    USE dener,                ONLY : denl, denl6, dekin6, detot
    USE uspp,                 ONLY : vkb, becsum, deeq, nkb, okvan, nlcc_any
    USE io_global,            ONLY : stdout, ionode
    USE core,                 ONLY : rhoc
    USE gvecw,                ONLY : ngw
    USE gvect,                ONLY : gg, g
    USE gvect,                ONLY : gstart, mill, eigts1, eigts2, eigts3
    USE cp_electronic_mass,   ONLY : emass
    USE efield_module,        ONLY : tefield, efield_berry_setup, berry_energy, &
                                     tefield2, efield_berry_setup2, berry_energy2
    USE cg_module,            ONLY : tcg
    USE ensemble_dft,         ONLY : tens, compute_entropy
    USE cp_interfaces,        ONLY : runcp_uspp, runcp_uspp_force_pairing, &
                                     strucf, phfacs, nlfh, vofrho, nlfl_bgrp, prefor
    USE cp_interfaces,        ONLY : rhoofr, ortho, wave_rand_init, elec_fakekine
    USE cp_interfaces,        ONLY : compute_stress, dotcsc, calbec, caldbec_bgrp
    USE cp_interfaces,        ONLY : nlfq_bgrp
    USE printout_base,        ONLY : printout_pos
    USE orthogonalize_base,   ONLY : updatc, calphi_bgrp
    USE wave_base,            ONLY : wave_steepest
    USE wavefunctions,        ONLY : c0_bgrp, cm_bgrp, c0_d, phi, cm_d
    USE fft_base,             ONLY : dfftp, dffts
    USE time_step,            ONLY : delt
    USE cp_main_variables,    ONLY : idesc, bephi, becp_bgrp, nfi, iabox, nabox, &
                                     sfac, eigr, taub, irb, eigrb, bec_bgrp, bec_d, &
                                     lambda, lambdam, lambdap, ema0bg, rhog, rhor, rhos, &
                                     vpot, ht0, edft, becdr_bgrp, dbec, drhor, drhog
    USE mp_global,            ONLY : inter_bgrp_comm, nbgrp, me_bgrp
    USE mp_world,             ONLY : mpime, world_comm
    USE mp,                   ONLY : mp_sum, mp_barrier
    USE matrix_inversion
    USE device_memcpy_m,        ONLY : dev_memcpy

#if defined (__ENVIRON)
    USE plugin_flags,         ONLY : use_environ
    USE environ_base_module,  ONLY : update_environ_ions
#endif

    !
    IMPLICIT NONE
    !
    include 'laxlib.fh'
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
       CALL randpos( taus, nat, ityp, tranp, amprp, ainv, iforce )
       !
       CALL s_to_r( taus, tau0, nat, h )
       !
    END IF
    !
    CALL phfacs( eigts1, eigts2, eigts3, eigr, mill, taus, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat )
    !
    CALL strucf( sfac, eigts1, eigts2, eigts3, mill, dffts%ngm )
    !     
    IF ( okvan .OR. nlcc_any ) THEN
       CALL initbox ( tau0, alat, at, ainv, taub, irb, iabox, nabox )
       CALL phbox( taub, iverbosity, eigrb )
    END IF
    !
    !     pass ions informations to plugins
    !
#if defined(__LEGACY_PLUGINS)
  CALL plugin_init_ions(tau0)
#endif
#if defined (__ENVIRON)
    IF (use_environ) CALL update_environ_ions(tau0)
#endif
    !
    !     wfc initialization with random numbers
    !     
    CALL wave_rand_init( cm_bgrp )
    !
    IF ( ionode ) &
       WRITE( stdout, fmt = '(//,3X, "Wave Initialization: random initial wave-functions" )' )

    ! if asked, use as many atomic wavefunctions as possible
    if ( trim(startingwfc) == 'atomic') then
       if ( ionode ) &
         WRITE (stdout, '("Using also atomic wavefunctions as much as possible")') 
       CALL atomic_wfc_cp(omega, nat, nsp, ityp, tau0, iupdwn, nspin, &
               ngw, nbspx, cm_bgrp )
    endif


    !
    ! ... prefor calculates vkb (used by gram)
    !
    CALL prefor( eigr, vkb )
    !
    !$acc update device(vkb)
    !
    nspin_wfc = nspin
    IF( force_pairing ) nspin_wfc = 1

    CALL gram_bgrp( vkb, bec_bgrp, nkb, cm_bgrp, ngw )

    IF( force_pairing ) cm_bgrp(:,iupdwn(2):iupdwn(2)+nupdwn(2)-1) = cm_bgrp(:,1:nupdwn(2))
    !
    if( iverbosity > 1 ) CALL dotcsc( vkb, cm_bgrp, ngw, nbsp )
    !
#if defined(__CUDA)
    CALL dev_memcpy( cm_d, cm_bgrp )
#endif
    !
    ! ... initialize bands
    !
    CALL occn_info( f )
    !
    hold = h
    velh = 0.0d0
    fion = 0.0d0
    !
    IF ( tv0rd .AND. tfor .AND. .NOT. tcap ) THEN
       !
       ! ... vel_srt=starting velocities, read from input, are brough to
       ! ... scaled axis and copied into array vels. Since velocites are
       ! ... not actually used by the Verlet algorithm, we set tau(t-dt)
       ! ... to tausm=tau(t)-v*delta t so that the Verlet algorithm will 
       ! ... start with the correct velocity
       !

       CALL r_to_s( vel, vels, nat, ainv )
       tausm(:,:) =  taus(:,:) - vels(:,:)*delt
       velsm(:,:) =  vels(:,:)
    ELSE IF (tcap .AND. tfor ) THEN
       WRITE( stdout, '(" Randomizing ions velocities according to tempw (INPUT VELOCITIES DISCARDED) ")')
       CALL randvel( tempw, tau0 , taum, nat, ityp, iforce, amass, delt )
       CALL r_to_s( taum, tausm, nat, ainv )
       vels(:,:) = (taus(:,:)-tausm(:,:))/delt
       velsm(:,:) = vels(:,:)
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
       CALL calbec ( nbsp_bgrp, vkb, cm_bgrp, bec_bgrp, 0 )
       !
       if ( tstress ) CALL caldbec_bgrp( eigr, cm_bgrp, dbec, idesc )
       !
       CALL rhoofr( nfi, cm_bgrp, cm_d, bec_bgrp, dbec, becsum, rhor, drhor, rhog, drhog, rhos, enl, denl, ekin, dekin6 )
       !
       edft%enl  = enl
       edft%ekin = ekin
       !
    END IF
    !
    !     put core charge (if present) in rhoc(r)
    !
    if ( nlcc_any ) CALL set_cc( rhoc )
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
             CALL printout_pos( stdout, fion, nat, ityp, head = ' fion ' )

      CALL newd( vpot, becsum, fion, .true. )
      !
      IF( force_pairing ) THEN
         !
         CALL runcp_uspp_force_pairing( nfi, fccc, ccc, ema0bg, dt2bye, rhos,&
                    bec_bgrp, cm_bgrp, c0_bgrp, ei_unp, fromscra = .TRUE. )
         !
         CALL setval_lambda( lambda(:,:,2), nupdwn(1), nupdwn(1), 0.d0, idesc(:,1) )
         !
      ELSE
         !
         CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec_bgrp, cm_bgrp, cm_d, c0_bgrp, c0_d, fromscra = .TRUE. )
         !
      ENDIF
      !
#if defined(__CUDA)
      CALL dev_memcpy( c0_d, c0_bgrp )  ! c0 contains the updated wave functions
#endif
      !
      !     nlfq needs deeq bec
      !
      IF( ttforce ) THEN
#if defined (__CUDA)
         !$acc data present(vkb)
         !$acc host_data use_device(vkb)
         CALL nlfq_bgrp( cm_d, vkb, bec_bgrp, becdr_bgrp, fion )
         !$acc end host_data 
         !$acc end data 
#else
         CALL nlfq_bgrp( cm_bgrp, vkb, bec_bgrp, becdr_bgrp, fion )
#endif
      END IF
      !
      !     calphi calculates phi
      !     the electron mass rises with g**2
      !
#if defined (__CUDA)
      !$acc data present(vkb)
      !$acc host_data use_device(vkb)
      CALL calphi_bgrp( cm_d, ngw, bec_bgrp, nkb, vkb, phi, nbspx_bgrp, ema0bg )
      !$acc end host_data 
      !$acc end data 
#else
      CALL calphi_bgrp( cm_bgrp, ngw, bec_bgrp, nkb, vkb, phi, nbspx_bgrp, ema0bg )
#endif
      !
      IF( force_pairing ) THEN
         !   phi( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) =    phi( :, 1:nupdwn(2))
         CALL dev_memcpy(phi(:,iupdwn(2):), phi(:,1:),  [1, ngw], 1 , [1, nupdwn(2)], 1)
      END IF
      !
      if( tortho ) then
#if defined (__CUDA)
         !$acc data present(vkb)
         !$acc host_data use_device(vkb)
         CALL ortho( vkb, c0_d, phi, lambda, idesc, bigr, iter, ccc, bephi, becp_bgrp )
         !$acc end host_data 
         !$acc end data 
#else
         CALL ortho( vkb, c0_bgrp, phi, lambda, idesc, bigr, iter, ccc, bephi, becp_bgrp )
#endif
      else
         CALL gram_bgrp( vkb, bec_bgrp, nkb, c0_bgrp, ngw )
      endif
      !
      IF ( ttforce ) THEN
         CALL nlfl_bgrp( bec_bgrp, becdr_bgrp, lambda, idesc, fion )
      END IF

      if ( iverbosity > 1 ) &
         CALL laxlib_print_matrix( lambda, idesc, nbsp, 9, nudx, ccc, ionode, stdout )

      !
      if ( tstress ) CALL nlfh( stress, bec_bgrp, dbec, lambda, idesc )
      !
      IF ( tortho ) THEN
#if defined (__CUDA)
         CALL updatc( ccc, lambda, phi, bephi, becp_bgrp, bec_d, c0_d, idesc )
         CALL dev_memcpy( c0_bgrp, c0_d )
         CALL dev_memcpy( bec_bgrp, bec_d )
#else
         CALL updatc( ccc, lambda, phi, bephi, becp_bgrp, bec_bgrp, c0_bgrp, idesc )
#endif
      END IF
      !
      IF( force_pairing ) THEN
         !
         c0_bgrp ( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) = c0_bgrp( :, 1:nupdwn(2))
         CALL dev_memcpy(phi(:,iupdwn(2):), phi(:,1:),  [1, ngw], 1 , [1, nupdwn(2)], 1)
         !phi( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) = phi( :, 1:nupdwn(2))
         lambda(:,:,2) = lambda(:,:,1)
         !
      ENDIF
      !
      ! the following compute only on NC pseudo components
      CALL calbec ( nbsp_bgrp, vkb, c0_bgrp, bec_bgrp, 1 )
      !
      if ( tstress ) CALL caldbec_bgrp( eigr, cm_bgrp, dbec, idesc )

      if ( iverbosity > 1 ) CALL dotcsc( vkb, c0_bgrp, ngw, nbsp_bgrp )
      !
      xnhp0 = 0.0d0
      xnhpm = 0.0d0
      vnhp  = 0.0d0
      fionm = 0.0d0

      ! Is this call useless and wrong?
      CALL ions_vel( vels, taus, tausm, delt )
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

subroutine hangup
    USE mp_world,             ONLY : mpime, world_comm
    USE mp,                   ONLY : mp_sum, mp_barrier
    call mp_barrier(world_comm)
    CALL stop_cp_run()
end subroutine

SUBROUTINE atomic_wfc_cp(omega, nat, nsp, ityp, tau, iupdwn, npol, npw, nbspx,&
               evc )
   
         USE kinds,        ONLY : DP
         USE uspp_param,   ONLY : nwfcm
         USE mp_global,    ONLY : intra_bgrp_comm
         USE uspp_data,    ONLY : nqx, tab_at, dq
         USE gvecw,        ONLY : ecutwfc
         USE upf_ions,     ONLY : n_atom_wfc

         IMPLICIT NONE
         !
         INTEGER, INTENT(IN) :: nat, nsp, ityp(nat), iupdwn(2), npol, npw, nbspx
         REAL(DP), INTENT(IN) :: omega, tau(3,nat)
         COMPLEX(DP), INTENT(inout) :: evc (npw,nbspx)
         !
         INTEGER :: natomwfc
         COMPLEX(DP), ALLOCATABLE  :: wfcatom(:,:,:)
         !! Superposition of atomic wavefunctions
   
         !cp specific settings (gamma only)
         ! xk is 0,0,0, igk is the identical permutation
         ! wfcatom has 2 dimensions (npw, nbnd*nspin)
         ! looks like only nspin=1,2 are implemented. The layout of the wfc is different:
         ! in cp it is the equivalent of (npw, nbnd, nspin ), while in pw is (npwx, nspin, nbnd)
         real(dp) :: xk(3)
         integer ::  igk(npw), i, ipol, sh(2)
         REAL(DP) :: angle1(nsp), angle2(nsp)
         !! dummy, unused variables
   
         ! identity permutation
         do i=1,npw
            igk(i)=i
         enddo
         ! gamma point only
         xk=0.d0
         nqx = INT( (SQRT(ecutwfc) / dq + 4) * 1.d0 )
         allocate(tab_at(nqx,nwfcm,nsp))
         call init_tab_atwfc(omega, intra_bgrp_comm)

         natomwfc = n_atom_wfc ( nat, ityp )
         allocate ( wfcatom(npw, npol, natomwfc) )
         
         ! only nospin / LSDA in CP
         call atomic_wfc_acc( xk, npw, igk, nat, nsp, ityp, tau, &
              .false., .false.,  angle1, angle2, .false., &
              npw, npol, natomwfc, wfcatom )
   
         sh = shape(evc)
   
         !write the result in the correct order in evc
         do i=1,min(natomwfc,sh(2)/npol)
            do ipol = 1, npol
               evc(:,i + iupdwn(ipol)-1) = wfcatom(:,ipol,i)
            enddo
         enddo
   
         deallocate (wfcatom)
         deallocate (tab_at)
         
      end subroutine atomic_wfc_cp
