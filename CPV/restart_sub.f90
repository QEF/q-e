!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE from_restart_module
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: from_restart
  !
  INTERFACE from_restart
     !
     MODULE PROCEDURE from_restart_cp, &
                      from_restart_fpmd, &
                      from_restart_sm
     !
  END INTERFACE from_restart
  !
  CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE from_restart_cp( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst, &
                              fion, taub, irb, eigrb, b1, b2, b3, nfi, rhog, &
                              rhor, rhos, rhoc, stress, detot, enthal,       &
                              lambda, lambdam, lambdap, ema0bg, dbec, bephi, &
                              becp, velh, dt2bye, fionm, ekincm )
    !--------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE control_flags,        ONLY : tranp, trane, trhor, iprsta, tpre, &
                                     tzeroc, tzerop, tzeroe, tfor, thdyn, &
                                     lwf,  tprnfor, tortho, tv0rd, amprp, &
                                     taurdr, ortho_eps, ortho_max, nbeg, dt_old
    USE ions_positions,       ONLY : taus, tau0, tausm, taum, vels, velsm, &
                                     ions_hmove
    USE ions_base,            ONLY : na, nat, nsp, randpos, zv, ions_vel, &
                                     pmass, iforce, vel
    USE cell_base,            ONLY : ainv, h, s_to_r, ibrav, omega, press, &
                                     hold, r_to_s, deth, wmass, iforceh,   &
                                     cell_force, cell_hmove
    USE efcalc,               ONLY : ef_force
    USE electrons_base,       ONLY : nbsp, nspin, nupdwn, iupdwn
    USE uspp,                 ONLY : vkb, becsum, deeq, nkb
    USE wavefunctions_module, ONLY : c0, cm, phi => cp
    USE io_global,            ONLY : stdout
    USE cpr_subroutines,      ONLY : compute_stress
    USE wannier_subroutines,  ONLY : get_wannier_center, &
                                     write_charge_and_exit, &
                                     ef_tune, wf_options, ef_potential
    USE core,                 ONLY : nlcc_any
    USE gvecw,                ONLY : ngw
    USE gvecs,                ONLY : ngs
    USE reciprocal_vectors,   ONLY : gstart, mill_l
    USE wave_base,            ONLY : wave_verlet
    USE cvan,                 ONLY : nvb
    USE ions_nose,            ONLY : xnhp0, xnhpm, vnhp
    USE cp_electronic_mass,   ONLY : emass
    USE efield_module,        ONLY : efield_berry_setup, tefield, &
                                     efield_berry_setup2, tefield2
    USE runcp_module,         ONLY : runcp_uspp, runcp_uspp_force_pairing, runcp_uspp_bgl
    USE wave_constrains,      ONLY : interpolate_lambda
    USE energies,             ONLY : eself, enl, etot, ekin
    USE time_step,            ONLY : delt, tps
    USE electrons_nose,       ONLY : xnhe0, xnhem, vnhe
    USE cell_nose,            ONLY : xnhh0, xnhhm, vnhh, cell_nosezero
    USE phase_factors_module, ONLY : strucf
    USE cg_module,            ONLY : tcg
    USE orthogonalize,        ONLY : ortho
    USE orthogonalize_base,   ONLY : updatc, calphi
    use charge_density,       only : rhoofr
    USE control_flags,        ONLY : force_pairing
    USE wave_functions,       ONLY : elec_fakekine
    !
    COMPLEX(DP) :: eigr(:,:), ei1(:,:), ei2(:,:), ei3(:,:)
    COMPLEX(DP) :: eigrb(:,:)
    INTEGER     :: irb(:,:)
    REAL(DP)    :: bec(:,:), fion(:,:), becdr(:,:,:), fionm(:,:)
    REAL(DP)    :: taub(:,:)
    REAL(DP)    :: b1(:), b2(:), b3(:)
    INTEGER     :: nfi
    LOGICAL     :: tfirst
    COMPLEX(DP) :: sfac(:,:)
    COMPLEX(DP) :: rhog(:,:)
    REAL(DP)    :: rhor(:,:), rhos(:,:), rhoc(:)
    REAL(DP)    :: stress(:,:), detot(:,:), enthal, ekincm
    REAL(DP)    :: lambda(:,:,:), lambdam(:,:,:), lambdap(:,:,:)
    REAL(DP)    :: ema0bg(:)
    REAL(DP)    :: dbec(:,:,:,:)
    REAL(DP)    :: bephi(:,:), becp(:,:)
    REAL(DP)    :: velh(3,3)
    REAL(DP)    :: dt2bye
    !
    REAL(DP),    ALLOCATABLE :: emadt2(:), emaver(:)
    COMPLEX(DP), ALLOCATABLE :: c2(:), c3(:)
    REAL(DP)                 :: verl1, verl2
    REAL(DP)                 :: bigr
    INTEGER                  :: i, j, iter, iss
    LOGICAL                  :: tlast = .FALSE.
    REAL(DP)                 :: fcell(3,3)
    REAL(DP)                 :: fccc = 0.5D0
    REAL(DP)                 :: ccc  = 0.D0
    ! Kostya: the variable below will disable the ionic & cell motion
    ! which nobody has any clue about ...
    REAL(DP)                 :: delt0 = 0.D0
    REAL(DP)                 :: ei_unp 
    INTEGER                  :: n_spin_start 

    !
    !
    ! ... We are restarting from file recompute ainv
    !
    CALL invmat( 3, h, ainv, deth )
    !
    ! ... Reset total time counter if the run is not strictly 'restart'
    !
    IF ( nbeg < 1 ) tps = 0.D0
    !
    IF ( taurdr ) THEN
       !
       ! ... Input positions read from input file and stored in tau0
       !
       CALL r_to_s( tau0, taus, na, nsp, ainv )
       !
    END IF
    !
    IF ( ANY( tranp(1:nsp) ) ) THEN
       !
       ! ... Input positions are randomized
       !
       CALL randpos( taus, na, nsp, tranp, amprp, ainv, iforce )
       !
    END IF
    !
    CALL s_to_r( taus, tau0, na, nsp, h )
    !
    IF ( tzerop ) THEN
       !
       tausm = taus
       vels  = 0.D0
       !
    END IF
    !
    IF ( tzeroc ) THEN
       !
       hold = h
       velh = 0.D0
       !
    END IF
    !
    fion = 0.D0
    !
    IF ( tzeroe ) THEN
       !
       lambdam = lambda
       !
       IF( force_pairing ) c0(:,iupdwn(2):nbsp,1) = c0(:,1:nupdwn(2),1)
       !
       cm(:,1:nbsp,1) = c0(:,1:nbsp,1)
       !
       WRITE( stdout, '(" Electronic velocities set to zero")' )
       !
    END IF
    !
    CALL phfac( tau0, ei1, ei2, ei3, eigr )
    !
    CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
    !
    CALL prefor( eigr, vkb )
    !
    IF ( trane .AND. trhor ) THEN
       !
       CALL gram( vkb, bec, nkb, c0, ngw, nbsp )
       !
       IF( force_pairing ) c0(:,iupdwn(2):nbsp,1) = c0(:,1:nupdwn(2),1)
       !
       cm(:,1:nbsp,1) = c0(:,1:nbsp,1)
       !
    END IF
    !
    CALL formf( tfirst, eself )
    !
    CALL calbec( 1, nsp, eigr, c0, bec )
    !
    IF ( tpre ) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, c0, dbec, .true. )
    !
    IF ( tefield ) CALL efield_berry_setup( eigr, tau0 )
    IF ( tefield2 ) CALL efield_berry_setup2( eigr, tau0 )
    !
    IF ( tzerop .OR. tzeroe .OR. taurdr ) THEN
       !
       CALL initbox( tau0, taub, irb )
       !
       CALL phbox( taub, eigrb )
       !
       IF ( lwf ) &
          CALL get_wannier_center( tfirst, c0, bec, becdr, eigr, &
                                   eigrb, taub, irb, ibrav, b1, b2, b3 )
       !
       CALL rhoofr( nfi, c0(:,:,1), irb, eigrb, bec, &
                    becsum, rhor, rhog, rhos, enl, ekin )
       !
       ! ... put core charge (if present) in rhoc(r)
       !
       IF ( nlcc_any ) CALL set_cc( irb, eigrb, rhoc )
       !
       IF ( lwf ) THEN
          !
          CALL write_charge_and_exit( rhog )
          CALL ef_tune( rhog, tau0 )
          !
       END IF
       !
       IF ( .NOT. tcg ) THEN
          CALL vofrho( nfi, rhor, rhog, rhos, rhoc, tfirst, tlast, &
                       ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )
       !
       CALL compute_stress( stress, detot, h, omega )
       !
       IF ( lwf ) &
          CALL wf_options( tfirst, nfi, c0, becsum, bec, becdr, &
                           eigr, eigrb, taub, irb, ibrav, b1,   &
                           b2, b3, rhor, rhog, rhos, enl, ekin )
       !
       CALL newd( rhor, irb, eigrb, becsum, fion )
       !
       CALL prefor( eigr, vkb )
       
       IF ( tzeroe .AND. ( .NOT. tcg ) ) THEN 

         IF( force_pairing ) THEN
             !
             CALL runcp_uspp_force_pairing( nfi, fccc, ccc, ema0bg, dt2bye, rhos, &
                           bec, c0(:,:,1), cm(:,:,1), ei_unp, restart = .TRUE. )
             !          lambda(nupdwn(1), nupdwn(1), 1) = ei_unp
             !          lambda(nupdwn(1), nupdwn(1), 2) = 0.d0 
             lambdam( 1:nupdwn(2), 1:nupdwn(2), 2) = lambdam( 1:nupdwn(2), 1:nupdwn(2), 1)
             lambda( 1:nupdwn(2), 1:nupdwn(2), 2) =  lambda( 1:nupdwn(2), 1:nupdwn(2), 1)
             !
          ELSE
             !
#if defined __BGL
             CALL runcp_uspp_bgl( nfi, fccc, ccc, ema0bg, dt2bye, rhos, &
                           bec, c0(:,:,1), cm(:,:,1), restart = .TRUE. )
#else
             CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, &
                           bec, c0(:,:,1), cm(:,:,1), restart = .TRUE. )
#endif
             !
          ENDIF 

       ENDIF 


       !
       ! ... nlfq needs deeq bec
       !
       IF ( .NOT. tcg ) THEN
          !
          IF ( tfor .OR. tprnfor ) CALL nlfq( c0, eigr, bec, becdr, fion )
          !
          IF ( tfor .OR. thdyn ) then
             CALL interpolate_lambda( lambdap, lambda, lambdam )
          ELSE
             ! take care of the otherwise uninitialized lambdam
             lambdam = lambda
          END IF
          !
       END IF
       !
       ! ... calphi calculates phi; the electron mass rises with g**2
       !
       CALL calphi( c0, ngw, bec, nkb, vkb, phi, nbsp, ema0bg )
       !
       ! ... begin try and error loop ( only one step! )
       !
       ! ... nlfl and nlfh need: lambda (guessed) becdr
       !
       IF ( ( tfor .OR. tprnfor ) &
           .AND. .NOT. tcg ) CALL nlfl( bec, becdr, lambda, fion )
       !
       IF ( tpre ) CALL nlfh( bec, dbec, lambda )
       !
       IF ( tortho ) THEN
          !
          CALL ortho( eigr, cm(:,:,1), phi(:,:,1), ngw, lambda, SIZE(lambda,1), &
                      bigr, iter, dt2bye, bephi, becp, nbsp, nspin, nupdwn, iupdwn )
          !
          n_spin_start = nspin 
          IF( force_pairing ) n_spin_start = 1
          !
          DO iss = 1,n_spin_start !!nspin
             CALL updatc( dt2bye, nbsp, lambda(:,:,iss), SIZE(lambda,1), phi, SIZE(phi,1), &
                       bephi, SIZE(bephi,1), becp, bec, cm, nupdwn(iss),iupdwn(iss) )
          END DO
          !
       ELSE
          !
          IF( .not. tcg) CALL gram( vkb, bec, nkb, cm, ngw, nbsp )
          !
       END IF
       !
       IF( force_pairing ) THEN
            cm( :,iupdwn(2):(iupdwn(2)-1+nupdwn(2)),1) =     cm( :,1:nupdwn(2),1) 
           phi( :,iupdwn(2):(iupdwn(2)-1+nupdwn(2)),1) =    phi( :,1:nupdwn(2),1)
        lambda(1:nupdwn(2),1:nupdwn(2),2)                = lambda(1:nupdwn(2),1:nupdwn(2),1) 
       ENDIF 
       !
       CALL calbec( nvb+1, nsp, eigr, cm, bec )
       !
       IF ( tpre ) &
          CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, cm, dbec, .TRUE. )
       !
       IF ( thdyn ) THEN
          !
          CALL cell_force( fcell, ainv, stress, omega, press, wmass )
          !
          CALL cell_hmove( h, hold, delt0, iforceh, fcell )
          !
          CALL invmat( 3, h, ainv, deth )
          !
       END IF
       !
       IF ( tfor ) THEN
          !
          IF ( lwf ) CALL ef_force( fion, na, nsp, zv )
          !
          IF ( tv0rd ) THEN
             !
             CALL r_to_s( vel, vels, na, nsp, h )
             !
             taus(:,:) = tausm(:,:) + delt * vels(:,:)
             !
          ELSE
             !
             CALL ions_hmove( taus, tausm, iforce, &
                              pmass, fion, ainv, delt0, na, nsp )
             !
          END IF
          !
          CALL s_to_r( taus, tau0, na, nsp, h )
          !
          CALL phfac( tau0, ei1, ei2, ei3, eigr )
          !
          CALL calbec( 1, nsp, eigr, c0, bec )
          !
          IF ( tpre ) &
             CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, c0, dbec, .TRUE. )
          !
       END IF
       !
       xnhp0 = 0.D0
       xnhpm = 0.D0
       vnhp  = 0.D0
       fionm = 0.D0
       !
       CALL ions_vel( vels, taus, tausm, na, nsp, delt )
       !
       CALL s_to_r( vels, vel, na, nsp, h )
       !
       CALL cell_nosezero( vnhh, xnhh0, xnhhm )
       !
       velh = ( h - hold ) / delt
       !
       !
       !     kinetic energy of the electrons
       !
       !
       IF( force_pairing ) THEN
          cm( :,iupdwn(2):(iupdwn(2)-1+nupdwn(2)),1) =     cm( :,1:nupdwn(2),1) 
          c0( :,iupdwn(2):(iupdwn(2)-1+nupdwn(2)),1) =     c0( :,1:nupdwn(2),1) 
          lambda(1:nupdwn(2),1:nupdwn(2),2)                = lambda(1:nupdwn(2),1:nupdwn(2),1) 
       ENDIF 
       !
       lambdam = lambda
       !
       CALL elec_fakekine( ekincm, ema0bg, emass, c0, cm, ngw, nbsp, 1, delt )
       !
       xnhe0 = 0.D0
       xnhem = 0.D0
       vnhe  = 0.D0
       !
       CALL DSWAP( 2*ngw*nbsp, c0, 1, cm, 1 )
       !
       END IF
    END IF
    !
    ! dt_old should be -1.0 here if untouched ...
    !
    if (dt_old.gt.0.0d0) then
       tausm = taus - (taus-tausm)*delt/dt_old
       xnhpm = xnhp0 - (xnhp0-xnhpm)*delt/dt_old
       WRITE( stdout, '(" tausm & xnhpm were rescaled ")' )
    endif
    !
    RETURN
    !
  END SUBROUTINE from_restart_cp
  !
  !--------------------------------------------------------------------------
  SUBROUTINE from_restart_sm( tfirst, taus, tau0, h, eigr, &
                              bec, c0, cm, ei1, ei2, ei3, sfac, eself )
    !--------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE control_flags,        ONLY : trane, trhor, iprsta, tpre
    USE uspp,                 ONLY : vkb, nkb
    USE ions_base,            ONLY : na, nsp
    USE electrons_base,       ONLY : nbsp
    USE io_global,            ONLY : stdout
    USE cell_base,            ONLY : s_to_r
    USE cpr_subroutines,      ONLY : print_atomic_var
    USE reciprocal_vectors,   ONLY : gstart, mill_l
    USE gvecs,                ONLY : ngs
    USE gvecw,                ONLY : ngw
    USE phase_factors_module, ONLY : strucf
    USE cdvan,                ONLY : dbec
    !
    IMPLICIT NONE
    !
    LOGICAL     :: tfirst
    REAL(DP)    :: taus(:,:), tau0(:,:)
    REAL(DP)    :: h(3,3)
    COMPLEX(DP) :: eigr(:,:)
    REAL(DP)    :: bec(:,:)
    COMPLEX(DP) :: c0(:,:)
    COMPLEX(DP) :: cm(:,:)
    COMPLEX(DP) :: ei1(:,:)
    COMPLEX(DP) :: ei2(:,:)
    COMPLEX(DP) :: ei3(:,:)
    COMPLEX(DP) :: sfac(:,:)
    REAL(DP)    :: eself
    INTEGER     :: j
    !
    !
    CALL s_to_r( taus,  tau0, na, nsp, h )
    !
    IF ( iprsta > 2 ) THEN
       !
       CALL print_atomic_var( taus, na, nsp, ' read: taus ' )
       !
       WRITE( stdout, '(" read: cell parameters h ")' )
       WRITE( stdout, * )  ( h(1,j), j = 1, 3 )
       WRITE( stdout, * )  ( h(2,j), j = 1, 3 )
       WRITE( stdout, * )  ( h(3,j), j = 1, 3 )
       !
    END IF
    !
    CALL phfac( tau0, ei1, ei2, ei3, eigr )
    !
    CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
    !
    CALL prefor( eigr, vkb )
    !
    IF ( trane .AND. trhor ) THEN
       !
       CALL gram( vkb, bec, nkb, c0, ngw, nbsp )
       !
       cm(:,1:nbsp) = c0(:,1:nbsp)
       !
    END IF
    !
    CALL formf( tfirst, eself )
    !
    CALL calbec( 1, nsp, eigr, c0, bec )
    !
    IF ( tpre ) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, c0, dbec, .true. )
    !
    RETURN
    !
  END SUBROUTINE from_restart_sm
  !
  !--------------------------------------------------------------------------
  SUBROUTINE from_restart_fpmd( nfi, acc, rhoe, cm, c0, cdesc, &
                                eigr, ei1, ei2, ei3, sfac, fi, ht_m, ht_0, &
                                atoms_m, atoms_0, bec, becdr, vpot, edft )
    !--------------------------------------------------------------------------
    !
    ! ... this routine recreates the starting configuration from a 
    ! ... restart file
    !
    USE kinds,                 ONLY : DP
    USE phase_factors_module,  ONLY : strucf, phfacs
    USE time_step,             ONLY : delt
    USE charge_density,        ONLY : rhoofr
    USE wave_base,             ONLY : wave_verlet, rande_base
    USE electrons_module,      ONLY : occn_info
    USE electrons_base,        ONLY : nspin
    USE ions_base,             ONLY : na, nsp, nax, randpos, taui, cdmi
    USE ions_module,           ONLY : set_reference_positions, &
                                      print_scaled_positions, &
                                      set_velocities
    USE energies,              ONLY : dft_energy_type
    USE cell_module,           ONLY : boxdimensions, gethinv, alat
    USE cell_base,             ONLY : r_to_s, s_to_r
    USE nl,                    ONLY : nlrh_m
    USE potentials,            ONLY : vofrhos
    USE orthogonalize,         ONLY : ortho
    USE io_global,             ONLY : ionode, ionode_id
    USE io_global,             ONLY : stdout
    USE wave_types,            ONLY : wave_descriptor
    USE control_flags,         ONLY : tcarpar, nbeg, tranp, amprp, tfor, tsdp, &
                                      thdyn, tsdc, tbeg, tsde, tortho, tzeroe, &
                                      tzerop, tzeroc, taurdr, tv0rd, nv0rd,    &
                                      trane, ampre, force_pairing, iprsta,     &
                                      tprnfor, tpre
    USE parameters,            ONLY : nacx
    USE atoms_type_module,     ONLY : atoms_type
    USE ions_base,             ONLY : vel_srt, tau_units
    USE runcp_module,          ONLY : runcp_ncpp
    USE grid_dimensions,       ONLY : nr1, nr2, nr3
    USE reciprocal_vectors,    ONLY : mill_l
    USE gvecp,                 ONLY : ngm
    USE ions_base,             ONLY : nat, tau_srt
    USE uspp,                  ONLY : vkb, nkb
    !
    IMPLICIT NONE
    !
    INTEGER                    :: nfi
    REAL(DP)                   :: acc(nacx)
    COMPLEX(DP)                :: sfac(:,:)
    TYPE(atoms_type)           :: atoms_0, atoms_m
    COMPLEX(DP)                :: eigr(:,:)
    COMPLEX(DP)                :: ei1(:,:)
    COMPLEX(DP)                :: ei2(:,:)
    COMPLEX(DP)                :: ei3(:,:)
    COMPLEX(DP), INTENT(INOUT) :: cm(:,:,:), c0(:,:,:)
    REAL(DP)                   :: fi(:,:)
    TYPE(boxdimensions)        :: ht_m, ht_0
    REAL(DP)                   :: rhoe(:,:)
    TYPE(wave_descriptor)      :: cdesc
    REAL(DP)                   :: bec(:,:)
    REAL(DP)                   :: becdr(:,:,:)
    REAL(DP)                   :: vpot(:,:)
    TYPE(dft_energy_type)      :: edft
    !
    INTEGER     :: ig, ib, i, j, k, ik, nb, is, ia, ierr, isa, iss
    REAL(DP)    :: fccc
    REAL(DP)    :: stau(3), rtau(3), hinv(3,3)
    COMPLEX(DP) :: cgam(1,1,1)
    REAL(DP)    :: gam(1,1,1)
    LOGICAL     :: ttforce
    LOGICAL     :: tstress
    !
    !
    ! ... if tbeg .eq. true ht_0 is not read from the restart file, and
    ! ... has been already been initialized in subroutine init1 togheter
    ! ... with the g vectors modules
    ! ... if tbeg is false, ht_0 is read from the restart file and now
    ! ... we have to compute the inverse and the volume of the cell,
    ! ... together with the new reciprocal vectors
    !
    ttforce = tfor  .OR. tprnfor
    tstress = thdyn .OR. tpre
    fccc    = 0.5d0
    !
    IF ( .NOT. tbeg ) THEN
       !
       CALL newinit( ht_0%hmat )
       !
    END IF
    !
    CALL occn_info( fi )
    !
    IF ( taurdr ) THEN
       !
       ! ... positions are read from stdin (tau_srt) and not read from restart
       ! ... file then reinitialize structure atoms_0
       !
       hinv = TRANSPOSE( ht_0%m1 )
       !
       CALL r_to_s( tau_srt, atoms_0%taus, atoms_0%na, atoms_0%nsp, hinv )
       !
       atoms_0%taur(1:3,1:nat) = tau_srt(1:3,1:nat)
       !
    END IF
    !
    ! ... diagnostics
    !
    IF ( ionode ) THEN
       !
       WRITE( stdout, 100 )
       !
       IF ( .NOT. tbeg ) THEN
          WRITE( stdout, 110 )
       ELSE
          WRITE( stdout, 120 )
       END IF
       !
       IF ( .NOT. taurdr ) THEN
          WRITE( stdout, 130 )
       ELSE
          WRITE( stdout, 140 )
       END IF
       !
       IF ( tfor .AND. ( .NOT. tsdp ) ) THEN
          !
          IF ( .NOT. tv0rd ) THEN
             !
             IF ( .NOT. tzerop ) THEN
                WRITE( stdout, 150 )
             ELSE
                WRITE( stdout, 155 )
             END IF
             !
          ELSE
             !
             WRITE( stdout, 160 )
             !
          END IF
          !
       END IF
       !
       IF ( iprsta > 1 ) &
          CALL print_scaled_positions( atoms_0, 6, 'from restart module' )
       !
       IF ( .NOT. tsde ) THEN
          !
          IF ( .NOT. tzeroe ) THEN
             WRITE( stdout, 170 )
          ELSE
             WRITE( stdout, 180 )
          END IF
          !
       END IF
       !
       WRITE( stdout, * )
       !
    END IF
    !
    IF ( tzeroc ) THEN
       !
       ht_m      = ht_0
       ht_0%hvel = 0.D0
       !
    END IF
    !
    IF ( ANY( tranp ) ) THEN
       !
       hinv = TRANSPOSE( ht_0%m1 )
       !
       CALL randpos( atoms_0%taus, atoms_0%na, &
                     atoms_0%nsp, tranp, amprp, hinv, atoms_0%mobile )
       !
       CALL s_to_r( atoms_0%taus, atoms_0%taur, &
                    atoms_0%na, atoms_0%nsp, ht_0%hmat )
       !
    END IF
    !
    IF ( tzerop .AND. tfor ) THEN
       !
       ! ... set initial velocities
       !
       CALL set_velocities( atoms_m, atoms_0, vel_srt, ht_0, delt )
       !
    END IF
    !
    ! ... computes form factors and initializes nl-pseudop. according
    ! ... to starting cell (from ndr or again fort.10)
    !
    CALL phfacs( ei1, ei2, ei3, eigr, mill_l, &
                 atoms_0%taus, nr1, nr2, nr3, atoms_0%nat )
    !
    CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngm )
    !
    CALL prefor( eigr, vkb )
    !
    CALL formf( .TRUE. , edft%eself )
    !
    IF ( trane ) THEN
       !
       WRITE( stdout, 515 ) ampre
       !
       DO iss = 1, cdesc%nspin
          call rande_base( c0( :, :, iss), ampre )
          CALL gram( vkb, bec, nkb, c0(1,1,iss), SIZE(c0,1), cdesc%nbt( iss ) )
          call rande_base( cm( :, :, iss), ampre )
          CALL gram( vkb, bec, nkb, cm(1,1,iss), SIZE(cm,1), cdesc%nbt( iss ) )
       END DO
       !
    END IF
    !
    IF ( tzeroe .OR. tzerop ) THEN
       !
       ! ... set velocities to zero
       ! ... set right initial conditions when c0=cm or stau0=staum
       ! ... (the cell is kept fixed)
       !
       atoms_0%for = 0.D0
       ! 
       edft%enl = nlrh_m( c0, cdesc, ttforce, atoms_0%for, bec, becdr, eigr )
       !
       CALL rhoofr( nfi, c0, cdesc, fi, rhoe, ht_0 )
       !
       CALL vofrhos( .true. , ttforce, tstress, rhoe, &
                     atoms_0, vpot, bec, c0, cdesc, fi, eigr, &
                     ei1, ei2, ei3, sfac, ht_0, edft )
       !
       IF ( tzeroe ) THEN
          !
          IF ( tcarpar .AND. ( .NOT. force_pairing ) ) THEN
             !
             CALL runcp_ncpp( cm, cm, c0, cdesc, vpot, vkb, fi, bec, fccc, gam, cgam, restart = .TRUE. )
             !
             IF ( tortho ) THEN
                !
                CALL ortho( cm, c0, cdesc )
                !
             ELSE
                !
                DO iss = 1, nspin
                   ! 
                   CALL gram( vkb, bec, nkb, c0(1,1,iss), SIZE(c0,1), cdesc%nbt( iss ) )
                   ! 
                END DO
                !
             END IF
             !
          ELSE
             !
             cm = c0
             !
          END IF
          !
       END IF
       !
    END IF
    !
    ! ... reset some variables if nbeg < 0 
    ! ... ( new simulation or step counter reset to 0 )
    !
    IF ( nbeg <= 0 ) THEN
       !
       acc = 0.D0
       nfi = 0
       !
       CALL set_reference_positions( cdmi, taui, atoms_0, ht_0 )
       !
    END IF
    !  
    RETURN 
    ! 
100 FORMAT( /,3X,'MD PARAMETERS READ FROM RESTART FILE',/ &
             ,3X,'------------------------------------' )
110 FORMAT(   3X,'Cell variables From RESTART file' )
120 FORMAT(   3X,'Cell variables From INPUT file' )
130 FORMAT(   3X,'Ions positions From RESTART file' )
140 FORMAT(   3X,'Ions positions From INPUT file' )
150 FORMAT(   3X,'Ions Velocities From RESTART file' )
155 FORMAT(   3X,'Ions Velocities set to ZERO' )
160 FORMAT(   3X,'Ions Velocities From STANDARD INPUT' )
170 FORMAT(   3X,'Electronic Velocities From RESTART file' )
180 FORMAT(   3X,'Electronic Velocities set to ZERO' )
515 FORMAT(   3X,'Initial random displacement of el. coordinates',/ &
              3X,'Amplitude = ',F10.6 )
    !
  END SUBROUTINE from_restart_fpmd
  !
END MODULE from_restart_module
