!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!=----------------------------------------------------------------------------=!
MODULE from_scratch_module
!=----------------------------------------------------------------------------=!
  !
  IMPLICIT NONE
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: from_scratch
  !
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE from_scratch( sfac, eigr, ei1, ei2, ei3, bec, becdr,  &
                           taub, irb, eigrb, b1, b2, b3, nfi,  &
                           rhog, rhor, rhos, rhoc, lambda, lambdam, lambdap, ema0bg,  &
                           dbec, bephi, becp, atoms, edft, ht, vpot )
    !--------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE control_flags,        ONLY : tranp, trane, iprsta, tpre, tcarpar,  &
                                     tzeroc, tzerop, tzeroe, tfor, thdyn, &
                                     lwf, tprnfor, tortho, amprp, ampre,  &
                                     tsde, ortho_eps, ortho_max, program_name, &
                                     force_pairing, use_task_groups
    USE ions_positions,       ONLY : taus, tau0, tausm, vels, fion, fionm
    USE ions_base,            ONLY : na, nsp, randpos, zv, ions_vel, pmass
    USE ions_base,            ONLY : taui, cdmi, nat, iforce
    USE ions_nose,            ONLY : xnhp0, xnhpm, vnhp
    USE cell_base,            ONLY : ainv, h, s_to_r, ibrav, omega, press, &
                                     hold, r_to_s, deth, wmass, iforceh,   &
                                     cell_force, boxdimensions, velh
    USE cell_nose,            ONLY : xnhh0, xnhhm, vnhh
    USE electrons_nose,       ONLY : xnhe0, xnhem, vnhe
    use electrons_base,       ONLY : nbsp, f, nspin, nupdwn, iupdwn
    USE electrons_module,     ONLY : occn_info
    USE energies,             ONLY : entropy, eself, enl, ekin, enthal, etot, ekincm
    USE energies,             ONLY : dft_energy_type, debug_energies
    USE dener,                ONLY : denl, denl6, dekin6, detot
    USE uspp,                 ONLY : vkb, becsum, deeq, nkb
    USE io_global,            ONLY : stdout, ionode
    USE core,                 ONLY : nlcc_any
    USE gvecw,                ONLY : ngw
    USE gvecs,                ONLY : ngs
    USE gvecp,                ONLY : ngm
    USE reciprocal_vectors,   ONLY : gstart, mill_l, gx
    USE cvan,                 ONLY : nvb
    USE cp_electronic_mass,   ONLY : emass
    USE efield_module,        ONLY : tefield, efield_berry_setup, berry_energy, &
                                     tefield2, efield_berry_setup2, berry_energy2
    USE cg_module,            ONLY : tcg
    USE ensemble_dft,         ONLY : tens, compute_entropy
    USE cp_interfaces,        ONLY : runcp_uspp, runcp_uspp_force_pairing, &
                                     strucf, phfacs, nlfh
    USE cp_interfaces,        ONLY : rhoofr, ortho, nlrh, wave_rand_init, elec_fakekine
    USE cp_interfaces,        ONLY : vofrhos, compute_stress
    USE cp_interfaces,        ONLY : printout, print_lambda
    USE printout_base,        ONLY : printout_pos
    USE orthogonalize_base,   ONLY : updatc, calphi
    USE atoms_type_module,    ONLY : atoms_type
    USE wave_base,            ONLY : wave_steepest
    USE wavefunctions_module, ONLY : c0, cm, phi => cp
    USE grid_dimensions,      ONLY : nr1, nr2, nr3
    USE time_step,            ONLY : delt
    USE cp_main_variables,    ONLY : setval_lambda, descla
    USE mp_global,            ONLY : np_ortho, me_ortho, ortho_comm


    !
    IMPLICIT NONE
    !
    COMPLEX(DP) :: eigr(:,:), ei1(:,:),  ei2(:,:),  ei3(:,:)
    COMPLEX(DP) :: eigrb(:,:)
    REAL(DP)    :: bec(:,:), becdr(:,:,:)
    REAL(DP)    :: taub(:,:)
    REAL(DP)    :: b1(:), b2(:), b3(:)
    INTEGER     :: irb(:,:)
    INTEGER     :: nfi
    COMPLEX(DP) :: sfac(:,:)
    COMPLEX(DP) :: rhog(:,:)
    REAL(DP)    :: rhor(:,:), rhos(:,:), rhoc(:)
    REAL(DP)    :: lambda(:,:,:), lambdam(:,:,:), lambdap(:,:,:)
    REAL(DP)    :: ema0bg(:)
    REAL(DP)    :: dbec(:,:,:,:)
    REAL(DP)    :: bephi(:,:), becp(:,:)
    TYPE(atoms_type)    ,    INTENT(OUT)   :: atoms
    TYPE(dft_energy_type) ,  INTENT(OUT)   :: edft
    TYPE(boxdimensions)  ,   INTENT(INOUT) :: ht
    REAL(DP),                INTENT(OUT)   :: vpot(:,:)

    !
    REAL(DP),    ALLOCATABLE :: emadt2(:), emaver(:)
    COMPLEX(DP), ALLOCATABLE :: c2(:), c3(:)
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
       atoms%taus(:,1:nat) = taus(:,1:nat)
       atoms%taur(:,1:nat) = tau0(:,1:nat)
       !
    END IF
    !
    CALL phfacs( ei1, ei2, ei3, eigr, mill_l, atoms%taus, nr1, nr2, nr3, atoms%nat )
    !
    CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
    !     
    CALL initbox ( tau0, taub, irb )
    !     
    CALL phbox( taub, eigrb )
    !
    !     random initialization
    !     
    CALL wave_rand_init( cm, nbsp, 1 )
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

    DO iss = 1, nspin_wfc
       !
       CALL gram( vkb, bec, nkb, cm(1,iupdwn(iss)), ngw, nupdwn(iss) )
       !
       !DO i = iupdwn(iss), nupdwn( iss )
       !   DO j = iupdwn(iss), nupdwn( iss )
       !      CALL dotcsv( dum, eigr, cm(1,i), cm(1,j), ngw )
       !      WRITE(stdout,*) i, j, dum
       !   END DO
       !END DO
       !
    END DO

    IF( force_pairing ) cm(:,iupdwn(2):iupdwn(2)+nupdwn(2)-1) = cm(:,1:nupdwn(2))
    !
    if( iprsta >= 3 ) CALL dotcsc( eigr, cm, ngw, nbsp )
    !
    ! ... initialize bands
    !
    CALL occn_info( f )
    !
    atoms%for  = 0.D0
    atoms%vels = 0.D0
    !
    ! ... compute local form factors
    !
    hold = h
    velh = 0.0d0
    fion = 0.0d0
    tausm = taus

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
    !
    IF( program_name == 'CP90' ) THEN

       IF( .NOT. tcg ) THEN

         CALL calbec ( 1, nsp, eigr, cm, bec )
         if (tpre) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, cm, dbec, .true. )

         CALL initbox ( tau0, taub, irb )
         CALL phbox( taub, eigrb )
         !
         CALL rhoofr ( nfi, cm(:,:), irb, eigrb, bec, becsum, rhor, rhog, rhos, enl, denl, ekin, dekin6 )
         !
         !     put core charge (if present) in rhoc(r)
         !
         if ( nlcc_any ) CALL set_cc( irb, eigrb, rhoc )
   
         IF( tens ) THEN
           CALL compute_entropy( entropy, f(1), nspin )
           entropy = entropy * nbsp
         END IF
         !
         vpot = rhor
         !
         CALL vofrho( nfi, vpot, rhog, rhos, rhoc, tfirst, tlast, &
        &  ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )

         IF( tefield ) THEN
           CALL berry_energy( enb, enbi, bec, cm(:,:), fion ) 
           etot = etot + enb + enbi
         END IF
         IF( tefield2 ) THEN
           CALL berry_energy2( enb, enbi, bec, cm(:,:), fion )
           etot = etot + enb + enbi
         END IF

         CALL compute_stress( stress, detot, h, omega )

         if(iprsta.gt.2) &
             CALL printout_pos( stdout, fion, nat, head = ' fion ' )

         CALL newd( vpot, irb, eigrb, becsum, fion )
         !
         IF( force_pairing ) THEN
            !
            CALL runcp_uspp_force_pairing( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, cm, &
        &                 c0, ei_unp, fromscra = .TRUE. )
            !
            CALL setval_lambda( lambda(:,:,2), nupdwn(1), nupdwn(1), 0.d0, descla(:,1) )
            !
         ELSE
            !
            CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, cm, c0, fromscra = .TRUE. )
            !
         ENDIF
         !
         !     nlfq needs deeq bec
         !
         if( tfor .or. tprnfor ) CALL nlfq( cm, eigr, bec, becdr, fion )
         !
         !     calphi calculates phi
         !     the electron mass rises with g**2
         !
         CALL calphi( cm, ngw, bec, nkb, vkb, phi, nbsp, ema0bg )
         !
         IF( force_pairing ) &
         &   phi( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) =    phi( :, 1:nupdwn(2))


         if( tortho ) then
            CALL ortho( eigr, c0, phi, ngw, lambda, descla, &
                        bigr, iter, ccc, bephi, becp, nbsp, nspin, nupdwn, iupdwn )
         else
            CALL gram( vkb, bec, nkb, c0, ngw, nbsp )
         endif
         !
         !
         if ( tfor .or. tprnfor ) CALL nlfl( bec, becdr, lambda, fion )

         if ( iprsta >= 3 ) CALL print_lambda( lambda, nbsp, 9, ccc )

         if ( tpre ) CALL nlfh( stress, bec, dbec, lambda )
         !
         IF ( tortho ) THEN
            DO iss = 1, nspin_wfc
               CALL updatc( ccc, nbsp, lambda(:,:,iss), SIZE(lambda,1), phi, SIZE(phi,1), &
                            bephi, SIZE(bephi,1), becp, bec, c0, nupdwn(iss), iupdwn(iss), &
                            descla(:,iss) )
            END DO
         END IF
         !
         IF( force_pairing ) THEN
            !
            c0 ( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) =  c0( :, 1:nupdwn(2))
            phi( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) = phi( :, 1:nupdwn(2))
            lambda(:,:,2) = lambda(:,:,1)
            !
         ENDIF
         !
         CALL calbec ( nvb+1, nsp, eigr, c0, bec )

         if ( tpre ) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, cm, dbec, .true. )

         if ( iprsta >= 3 ) CALL dotcsc( eigr, c0, ngw, nbsp )
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
         CALL elec_fakekine( ekincm, ema0bg, emass, c0, cm, ngw, nbsp, 1, delt )

         xnhe0 = 0.0d0
         xnhem = 0.0d0
         vnhe  = 0.0d0

         lambdam = lambda
         !
       ELSE 
          !
          c0 = cm
          !
       END IF

    ELSE

       CALL nlrh( cm, ttforce, tstress, atoms%for, bec, becdr, eigr, edft%enl, denl6 )
       !
       CALL rhoofr( 0, tstress, cm, f, rhor, ht%deth, edft%ekin, dekin6 )
       !
       CALL vofrhos( ttprint, ttforce, tstress, rhor, atoms, &
                  vpot, bec, cm, f, eigr, ei1, ei2, ei3, sfac, ht, edft )
       !
       IF( iprsta > 1 ) CALL debug_energies( edft )
       !
       IF ( tcarpar ) THEN
          !
          IF ( .NOT. force_pairing ) THEN
             !
             CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, vpot, bec, cm, c0, fromscra = .TRUE. )
             !
          ELSE
             !
             c0 = cm
             !
          END IF
          !
          IF ( tortho .AND. ( .NOT. force_pairing ) ) THEN
             !
             ccc = fccc * dt2bye
             !
             CALL ortho( cm, c0, lambda, descla, ccc, nupdwn, iupdwn, nspin )
             !
          ELSE
             !
             DO iss = 1, nspin_wfc
               !
               CALL gram( vkb, bec, nkb, c0(1,iupdwn(iss)), SIZE(c0,1), nupdwn( iss ) )
               !
             END DO
             !
             IF( force_pairing ) c0(1,iupdwn(2):iupdwn(2)+nupdwn(2)-1) = c0(1,1:nupdwn(2))
             !
          END IF
          !
       ELSE 
          !
          c0 = cm
          !
       END IF
       !
    END IF

    IF( iprsta > 1 ) THEN
       !
       !  Printout values at step 0, useful for debugging
       !
       CALL printout( nfi, atoms, 0.0d0, 0.0d0, ttprint, ht, edft )
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE from_scratch
  !
!=----------------------------------------------------------------------------=!
END MODULE from_scratch_module
!=----------------------------------------------------------------------------=!
