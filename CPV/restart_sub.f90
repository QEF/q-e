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
                      from_restart_sm
     !
  END INTERFACE
  !
  CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE from_restart_cp( sfac, eigr, ei1, ei2, ei3, bec, becdr,         &
                              taub, irb, eigrb, b1, b2, b3, nfi, rhog,       &
                              rhor, rhos, rhoc,                              &
                              lambda, lambdam, lambdap, ema0bg, dbec, bephi, &
                              becp, htm, ht0, vpot, atoms0, edft )

    !--------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE control_flags,        ONLY : tranp, trane, trhor, iprsta, tpre, &
                                     tzeroc, tzerop, tzeroe, tfor, thdyn, &
                                     lwf,  tprnfor, tortho, tv0rd, amprp, &
                                     taurdr, ortho_eps, ortho_max, nbeg, dt_old, &
                                     use_task_groups, program_name, tsde, tcarpar
    USE ions_positions,       ONLY : taus, tau0, tausm, taum, vels, velsm, &
                                     ions_hmove, fion, fionm
    USE ions_base,            ONLY : na, nat, nsp, randpos, zv, ions_vel, &
                                     pmass, iforce, vel
    USE cell_base,            ONLY : ainv, h, s_to_r, ibrav, omega, press, &
                                     hold, r_to_s, deth, wmass, iforceh, velh,  &
                                     cell_force, cell_hmove, boxdimensions
    USE efcalc,               ONLY : ef_force
    USE electrons_base,       ONLY : nbsp, nspin, nupdwn, iupdwn, f
    USE uspp,                 ONLY : vkb, becsum, deeq, nkb
    USE wavefunctions_module, ONLY : c0, cm, phi => cp
    USE io_global,            ONLY : stdout
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
    USE cp_interfaces,        ONLY : runcp_uspp, runcp_uspp_force_pairing, &
                                     interpolate_lambda, from_restart_x, nlrh, vofrhos, &
                                     nlfh, phfacs
    USE energies,             ONLY : eself, enl, etot, ekin, dft_energy_type, enthal, ekincm
    USE time_step,            ONLY : delt, tps
    USE electrons_nose,       ONLY : xnhe0, xnhem, vnhe
    USE cell_nose,            ONLY : xnhh0, xnhhm, vnhh, cell_nosezero
    USE cg_module,            ONLY : tcg
    USE orthogonalize_base,   ONLY : updatc, calphi
    use cp_interfaces,        only : rhoofr, ortho, elec_fakekine, compute_stress
    USE control_flags,        ONLY : force_pairing
    USE dener,                ONLY : denl, dekin6, denl6, detot
    USE mp_global,            ONLY : np_ortho, me_ortho, ortho_comm, me_image
    USE cp_main_variables,    ONLY : descla
    USE atoms_type_module,    ONLY : atoms_type
    USE grid_dimensions,      ONLY : nr1, nr2, nr3
    !
    COMPLEX(DP) :: eigr(:,:), ei1(:,:), ei2(:,:), ei3(:,:)
    COMPLEX(DP) :: eigrb(:,:)
    INTEGER     :: irb(:,:)
    REAL(DP)    :: bec(:,:), becdr(:,:,:)
    REAL(DP)    :: taub(:,:)
    REAL(DP)    :: b1(:), b2(:), b3(:)
    INTEGER     :: nfi
    COMPLEX(DP) :: sfac(:,:)
    COMPLEX(DP) :: rhog(:,:)
    REAL(DP)    :: rhor(:,:), rhos(:,:), rhoc(:)
    REAL(DP)    :: lambda(:,:,:), lambdam(:,:,:), lambdap(:,:,:)
    REAL(DP)    :: ema0bg(:)
    REAL(DP)    :: dbec(:,:,:,:)
    REAL(DP)    :: bephi(:,:), becp(:,:)
    TYPE(boxdimensions)        :: htm, ht0
    REAL(DP)    :: vpot(:,:)
    TYPE(atoms_type) :: atoms0
    TYPE(dft_energy_type) :: edft
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
    REAL(DP)                 :: dt2bye 
    REAL(DP)                 :: stress(3,3) 
    INTEGER                  :: n_spin_start 
    LOGICAL                  :: ttforce
    LOGICAL                  :: tstress
    LOGICAL                  :: tfirst = .TRUE.
    !
    CALL from_restart_x( ht0, htm, phi, c0, cm, lambdap, lambda, lambdam, &
            ei1, ei2, ei3, eigr, sfac, vkb, nkb, bec, dbec, taub, irb, eigrb )
    !
    dt2bye = delt * delt / emass
    !
    stress = 0.0d0
    !
    IF( program_name == 'CP90' ) THEN
       !
       IF ( tzerop .OR. tzeroe .OR. taurdr ) THEN
          !
          IF ( lwf ) &
             CALL get_wannier_center( tfirst, c0, bec, becdr, eigr, &
                                      eigrb, taub, irb, ibrav, b1, b2, b3 )
          !
          CALL rhoofr( nfi, c0, irb, eigrb, bec, &
                       becsum, rhor, rhog, rhos, enl, denl, ekin, dekin6 )
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
          vpot = rhor
          !
          IF ( .NOT. tcg ) THEN
             CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
                          ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )
          !
          CALL compute_stress( stress, detot, h, omega )
          !
          IF ( lwf ) &
             CALL wf_options( tfirst, nfi, c0, becsum, bec, becdr, &
                              eigr, eigrb, taub, irb, ibrav, b1,   &
                              b2, b3, vpot, rhog, rhos, enl, ekin )
          !
          CALL newd( vpot, irb, eigrb, becsum, fion )
          !
          CALL prefor( eigr, vkb )
          
          IF ( tzeroe .AND. ( .NOT. tcg ) ) THEN 
   
            IF( force_pairing ) THEN
                !
                CALL runcp_uspp_force_pairing( nfi, fccc, ccc, ema0bg, dt2bye, rhos, &
                              bec, c0, cm, ei_unp, restart = .TRUE. )
                lambdam( :, :, 2) = lambdam( :, :, 1)
                lambda( :, :, 2) =  lambda( :, :, 1)
                !
             ELSE
                !
                CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, c0, cm, restart = .TRUE. )
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
          IF ( ( tfor .OR. tprnfor ) .AND. .NOT. tcg ) CALL nlfl( bec, becdr, lambda, fion )
          !
          IF ( tpre ) CALL nlfh( stress, bec, dbec, lambda )
          !
          IF ( tortho ) THEN
             !
             CALL ortho( eigr, cm, phi, ngw, lambda, descla, &
                         bigr, iter, dt2bye, bephi, becp, nbsp, nspin, nupdwn, iupdwn )
             !
             n_spin_start = nspin 
             IF( force_pairing ) n_spin_start = 1
             !
             DO iss = 1,n_spin_start !!nspin
                CALL updatc( dt2bye, nbsp, lambda(:,:,iss), SIZE(lambda,1), phi, SIZE(phi,1), &
                          bephi, SIZE(bephi,1), becp, bec, cm, nupdwn(iss),iupdwn(iss), &
                          descla(:,iss) )
             END DO
             !
          ELSE
             !
             IF( .not. tcg) CALL gram( vkb, bec, nkb, cm, ngw, nbsp )
             !
          END IF
          !
          IF( force_pairing ) THEN
               cm( :,iupdwn(2):(iupdwn(2)-1+nupdwn(2))) =     cm( :,1:nupdwn(2)) 
              phi( :,iupdwn(2):(iupdwn(2)-1+nupdwn(2))) =    phi( :,1:nupdwn(2))
           lambda(1:nupdwn(2),1:nupdwn(2),2)            = lambda(1:nupdwn(2),1:nupdwn(2),1) 
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
             CALL phfacs( ei1, ei2, ei3, eigr, mill_l, taus, nr1, nr2, nr3, nat )
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
             cm( :,iupdwn(2):(iupdwn(2)-1+nupdwn(2))) =     cm( :,1:nupdwn(2)) 
             c0( :,iupdwn(2):(iupdwn(2)-1+nupdwn(2))) =     c0( :,1:nupdwn(2)) 
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
          CALL DSWAP( 2*SIZE(c0), c0, 1, cm, 1 )
          !
          END IF
          !
       END IF
       !
       ! dt_old should be -1.0 here if untouched ...
       !
       if (dt_old.gt.0.0d0) then
          tausm = taus - (taus-tausm)*delt/dt_old
          xnhpm = xnhp0 - (xnhp0-xnhpm)*delt/dt_old
          WRITE( stdout, '(" tausm & xnhpm were rescaled ")' )
       endif

    ELSE
       !
       edft%eself = eself
       !
       IF ( taurdr .OR. ( tzeroe .AND. .NOT. tsde ) .OR. tzerop ) THEN
          !
          ttforce = tfor  .OR. tprnfor
          tstress = thdyn .OR. tpre
          fccc    = 0.5d0
          !
          ! ... set velocities to zero
          ! ... set right initial conditions when c0=cm or stau0=staum
          ! ... (the cell is kept fixed)
          !
          atoms0%for = 0.D0
          ! 
          CALL nlrh( c0, ttforce, tstress, atoms0%for, bec, becdr, eigr, edft%enl, denl6 )
          !
          CALL rhoofr( nfi, c0, irb, eigrb, bec, becsum, rhor, rhog, rhos, edft%enl, denl, edft%ekin, dekin6 )
          !
          CALL vofrhos( .true. , ttforce, tstress, rhor, rhog, &
                        atoms0, vpot, bec, c0, f, eigr, &
                        ei1, ei2, ei3, sfac, ht0, edft )
          !
          IF ( tzeroe ) THEN
             !
             IF ( tcarpar .AND. ( .NOT. force_pairing ) ) THEN
                !
                CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, vpot, bec, c0, cm, restart = .TRUE. )
                !
                !  now cm contains the unorthogonalized "c" at time t+dt
                !
                IF ( tortho ) THEN
                   !
                   ccc = fccc * dt2bye
                   !
                   CALL ortho( c0, cm, lambda, descla, ccc, nupdwn, iupdwn, nspin )
                   !
                ELSE
                   !
                   DO iss = 1, nspin
                      ! 
                      CALL gram( vkb, bec, nkb, cm(1,iupdwn(iss)), SIZE(cm,1), nupdwn( iss ) )
                      ! 
                   END DO
                   !
                END IF
                !
                CALL DSWAP( 2*SIZE(c0), c0, 1, cm, 1 ) 
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

       velh = htm%hvel
      
    END IF
    !
    RETURN
    !
  END SUBROUTINE from_restart_cp
  !
  !
  !--------------------------------------------------------------------------
  SUBROUTINE from_restart_sm( tfirst, taus, tau0, h, eigr, &
                              bec, c0, cm, ei1, ei2, ei3, sfac, eself )
    !--------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE control_flags,        ONLY : trane, trhor, iprsta, tpre
    USE uspp,                 ONLY : vkb, nkb
    USE ions_base,            ONLY : na, nsp, nat
    USE electrons_base,       ONLY : nbsp
    USE io_global,            ONLY : stdout
    USE cell_base,            ONLY : s_to_r
    USE printout_base,        ONLY : printout_pos
    USE reciprocal_vectors,   ONLY : gstart, mill_l
    USE gvecs,                ONLY : ngs
    USE gvecw,                ONLY : ngw
    USE cp_interfaces,        ONLY : strucf, phfacs
    USE cdvan,                ONLY : dbec
    USE grid_dimensions,      ONLY : nr1, nr2, nr3
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
       CALL printout_pos( stdout, taus, nat, head = ' read: taus ' )
       !
       WRITE( stdout, '(" read: cell parameters h ")' )
       WRITE( stdout, * )  ( h(1,j), j = 1, 3 )
       WRITE( stdout, * )  ( h(2,j), j = 1, 3 )
       WRITE( stdout, * )  ( h(3,j), j = 1, 3 )
       !
    END IF
    !
    CALL phfacs( ei1, ei2, ei3, eigr, mill_l, taus, nr1, nr2, nr3, nat )
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
  !
END MODULE from_restart_module




SUBROUTINE from_restart_x( &
   ht0, htm, phi, c0, cm, lambdap, lambda, lambdam, ei1, ei2, ei3, eigr, sfac, vkb, nkb, &
   bec, dbec, taub, irb, eigrb )
   !
   USE kinds,                 ONLY : DP
   USE control_flags,         ONLY : program_name, tbeg, taurdr, tfor, tsdp, tv0rd, &
                                     iprsta, tsde, tzeroe, tzerop, nbeg, tranp, amprp, &
                                     tzeroc, force_pairing, trhor, ampre, trane, tpre
   USE electrons_module,      ONLY : occn_info
   USE electrons_base,        ONLY : nspin, iupdwn, nupdwn, f, nbsp
   USE io_global,             ONLY : ionode, ionode_id, stdout
   USE cell_base,             ONLY : ainv, h, hold, deth, r_to_s, s_to_r, boxdimensions, &
                                     velh, a1, a2, a3
   USE ions_base,             ONLY : na, nsp, iforce, vel_srt, nat, randpos
   USE time_step,             ONLY : tps, delt
   USE ions_positions,        ONLY : taus, tau0, tausm, taum, vels, fion, fionm
   USE grid_dimensions,       ONLY : nr1, nr2, nr3
   USE reciprocal_vectors,    ONLY : mill_l
   USE printout_base,         ONLY : printout_pos
   USE gvecs,                 ONLY : ngs
   USE gvecw,                 ONLY : ngw
   USE cp_interfaces,         ONLY : phfacs, strucf
   USE energies,              ONLY : eself
   USE wave_base,             ONLY : rande_base
   USE mp_global,             ONLY : me_image
   USE efield_module,         ONLY : efield_berry_setup,  tefield, &
                                     efield_berry_setup2, tefield2
   USE small_box,             ONLY : ainvb
   USE uspp,                  ONLY : okvan
   !
   IMPLICIT NONE
   !
   TYPE (boxdimensions) :: ht0, htm
   COMPLEX(DP)          :: phi(:,:)
   COMPLEX(DP)          :: c0(:,:)
   COMPLEX(DP)          :: cm(:,:)
   REAL(DP)             :: lambda(:,:,:), lambdam(:,:,:), lambdap(:,:,:)
   COMPLEX(DP)          :: eigr(:,:)
   COMPLEX(DP)          :: ei1(:,:)
   COMPLEX(DP)          :: ei2(:,:)
   COMPLEX(DP)          :: ei3(:,:)
   COMPLEX(DP)          :: sfac(:,:)
   COMPLEX(DP)          :: vkb(:,:)
   INTEGER, INTENT(IN)  :: nkb
   REAL(DP)             :: bec(:,:)
   REAL(DP)             :: dbec(:,:,:,:)
   REAL(DP)             :: taub(:,:)
   INTEGER              :: irb(:,:)
   COMPLEX(DP)          :: eigrb(:,:)
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
      ! ... in readfile, only scaled positions are read
      ! ... file then reinitialize structure atoms_0
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
   IF ( tzerop .AND. tfor ) THEN
      !
      CALL r_to_s( vel_srt, vels, na, nsp, ainv )
      !
      CALL set_velocities( tausm, taus, vels, iforce, nat, delt )
      !
   END IF
   !
   CALL s_to_r( taus,  tau0, na, nsp, h )
   !
   CALL s_to_r( tausm, taum, na, nsp, h )
   !
   IF ( tzeroc ) THEN
      !
      hold = h
      velh = 0.D0
      !
      htm      = ht0
      ht0%hvel = 0.D0
      !
   END IF
   !
   fion = 0.D0
   !
   IF ( tzeroe ) THEN
      !
      IF( program_name == 'CP90' ) lambdam = lambda
      !
      IF( force_pairing ) c0(:,iupdwn(2):nbsp) = c0(:,1:nupdwn(2))
      !
      cm = c0
      !
      WRITE( stdout, '(" Electronic velocities set to zero")' )
      !
   END IF
   !
   IF( program_name == "FPMD" ) THEN
      !
      CALL occn_info( f )
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
            CALL printout_pos( stdout, taus, nat, head = 'Scaled positions from restart module' )
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
   END IF
   !
   ! ... computes form factors and initializes nl-pseudop. according
   ! ... to starting cell (from ndr or again standard input)
   !
   IF ( okvan) THEN
      CALL initbox( tau0, taub, irb, ainv, a1, a2, a3 )
      CALL phbox( taub, eigrb, ainvb )
   END IF
   !
   CALL phfacs( ei1, ei2, ei3, eigr, mill_l, taus, nr1, nr2, nr3, nat )
   !
   CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
   !
   CALL prefor( eigr, vkb )
   !
   CALL formf( .TRUE. , eself )
   !
   IF ( trane ) THEN
      !
      WRITE( stdout, 515 ) ampre
      !
515   FORMAT(   3X,'Initial random displacement of el. coordinates',/ &
                3X,'Amplitude = ',F10.6 )
      !
      CALL rande_base( c0, ampre )

      CALL gram( vkb, bec, nkb, c0, ngw, nbsp )
      !
      IF( force_pairing ) c0(:,iupdwn(2):nbsp) = c0(:,1:nupdwn(2))
      !
      cm = c0
      !
   END IF
   !

   IF( program_name == "CP90" ) THEN
      !
      CALL calbec( 1, nsp, eigr, c0, bec )
      !
      IF ( tpre ) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, c0, dbec, .true. )
      !
      IF ( tefield ) CALL efield_berry_setup( eigr, tau0 )
      IF ( tefield2 ) CALL efield_berry_setup2( eigr, tau0 )
      !
   END IF

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
    !
   !
   RETURN
END SUBROUTINE
