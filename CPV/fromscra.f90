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
  SUBROUTINE from_scratch( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst,   &
                              eself, fion, taub, irb, eigrb, b1, b2, b3, nfi,  &
                              rhog, rhor, rhos, rhoc, enl, ekin, stress, detot,&
                              enthal, etot, lambda, lambdam, lambdap, ema0bg,  &
                              dbec, delt, bephi, becp, velh, dt2bye, iforce,   &
                              fionm, xnhe0, xnhem, vnhe, ekincm, atoms, edft,  &
                              ht, cdesc, vpot )
    !--------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE parameters,           ONLY : nacx
    USE control_flags,        ONLY : tranp, trane, iprsta, tpre, tcarpar,  &
                                     tzeroc, tzerop, tzeroe, tfor, thdyn, &
                                     lwf, tprnfor, tortho, amprp, ampre,  &
                                     tsde, ortho_eps, ortho_max, program_name, &
                                     force_pairing
    USE ions_positions,       ONLY : taus, tau0, tausm, vels
    USE ions_base,            ONLY : na, nsp, randpos, zv, ions_vel, pmass
    USE ions_base,            ONLY : taui, cdmi, nat
    USE cell_base,            ONLY : ainv, h, s_to_r, ibrav, omega, press, &
                                     hold, r_to_s, deth, wmass, iforceh,   &
                                     cell_force
    USE cell_nose,            ONLY : xnhh0, xnhhm,  vnhh
    USE cell_module,          ONLY : boxdimensions
    use electrons_base,       ONLY : nbsp
    USE electrons_base,       ONLY : f, nspin, nupdwn, iupdwn
    USE electrons_module,     ONLY : occn_info
    USE energies,             ONLY : entropy
    USE energies,             ONLY : dft_energy_type, debug_energies
    USE uspp,                 ONLY : vkb, becsum, deeq, nkb
    USE io_global,            ONLY : stdout, ionode
    USE cpr_subroutines,      ONLY : compute_stress, print_atomic_var, &
                                     print_lambda
    USE core,                 ONLY : nlcc_any
    USE gvecw,                ONLY : ngw
    USE reciprocal_vectors,   ONLY : gstart, mill_l, gx
    USE gvecs,                ONLY : ngs
    USE gvecp,                ONLY : ngm
    USE cvan,                 ONLY : nvb
    USE ions_nose,            ONLY : xnhp0,  xnhpm,  vnhp
    USE cp_electronic_mass,   ONLY : emass
    USE efield_module,        ONLY : tefield, efield_berry_setup, &
                                     berry_energy, &
                                     tefield2, efield_berry_setup2, &
                                     berry_energy2
    USE cg_module,            ONLY : tcg
    USE ensemble_dft,         ONLY : tens, compute_entropy
    USE cp_interfaces,        ONLY : runcp_uspp, runcp_ncpp, runcp_uspp_force_pairing, &
                                     runcp_uspp_bgl, strucf, phfacs
    USE cp_interfaces,        ONLY : rhoofr, ortho, nlrh, wave_rand_init, elec_fakekine
    USE orthogonalize_base,   ONLY : updatc, calphi
    USE atoms_type_module,    ONLY : atoms_type
    USE wave_base,            ONLY : wave_steepest
    USE wavefunctions_module, ONLY : c0, cm, phi => cp
    USE wave_types,           ONLY : wave_descriptor
    USE potentials,           ONLY : vofrhos
    USE grid_dimensions,      ONLY : nr1, nr2, nr3
    USE cp_interfaces,        ONLY : printout


    !
    IMPLICIT NONE
    !
    COMPLEX(DP) :: eigr(:,:), ei1(:,:),  ei2(:,:),  ei3(:,:)
    COMPLEX(DP) :: eigrb(:,:)
    REAL(DP)    :: bec(:,:), fion(:,:), becdr(:,:,:), fionm(:,:)
    REAL(DP)    :: eself
    REAL(DP)    :: taub(:,:)
    REAL(DP)    :: b1(:), b2(:), b3(:)
    INTEGER     :: irb(:,:)
    INTEGER     :: nfi, iforce(:,:)
    LOGICAL     :: tfirst
    COMPLEX(DP) :: sfac(:,:)
    COMPLEX(DP) :: rhog(:,:)
    REAL(DP)    :: rhor(:,:), rhos(:,:), rhoc(:), enl, ekin
    REAL(DP)    :: stress(:,:), detot(:,:), enthal, etot
    REAL(DP)    :: lambda(:,:,:), lambdam(:,:,:), lambdap(:,:,:)
    REAL(DP)    :: ema0bg(:)
    REAL(DP)    :: dbec(:,:,:,:)
    REAL(DP)    :: delt
    REAL(DP)    :: bephi(:,:), becp(:,:)
    REAL(DP)    :: velh(:,:)
    REAL(DP)    :: dt2bye, xnhe0, xnhem, vnhe, ekincm
    TYPE(atoms_type)    ,    INTENT(OUT)   :: atoms
    TYPE(dft_energy_type) ,  INTENT(OUT)   :: edft
    TYPE(boxdimensions)  ,   INTENT(INOUT) :: ht
    TYPE(wave_descriptor),   INTENT(IN)    :: cdesc
    REAL(DP),                INTENT(OUT)   :: vpot(:,:)

    !
    REAL(DP),    ALLOCATABLE :: emadt2(:), emaver(:)
    COMPLEX(DP), ALLOCATABLE :: c2(:), c3(:)
    REAL(DP)                 :: verl1, verl2
    REAL(DP)                 :: bigr, dum
    REAL(DP)                 :: adum( nacx )
    INTEGER                  :: i, j, iter, iss, ierr, nspin_wfc
    LOGICAL                  :: tlast = .FALSE.
    REAL(DP)                 :: gam(1,1,1)
    REAL(DP)                 :: fcell(3,3), ccc, enb, enbi, fccc
    LOGICAL                  :: ttforce
    LOGICAL                  :: tstress
    LOGICAL, PARAMETER       :: ttprint = .TRUE.
    REAL(DP)                 :: ei_unp  
    INTEGER                  :: n_spin_start 
    !
    ttforce = tfor  .or. tprnfor
    tstress = thdyn .or. tpre
    !
    IF( tsde ) THEN
       fccc = 1.0d0
    ELSE
       fccc = 0.5d0
    END IF
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
    END DO

    IF( force_pairing ) cm(:,iupdwn(2):iupdwn(2)+nupdwn(2)-1) = cm(:,1:nupdwn(2))
    !
    if( iprsta .ge. 3 ) CALL dotcsc( eigr, cm, ngw, nbsp )
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

    IF( program_name == 'CP90' ) THEN
       lambdam = lambda
    END IF
    
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
         CALL rhoofr ( nfi, cm(:,:), irb, eigrb, bec, becsum, rhor, rhog, rhos, enl, ekin )
         !
         !     put core charge (if present) in rhoc(r)
         !
         if ( nlcc_any ) CALL set_cc( irb, eigrb, rhoc )
   
         IF( tens ) THEN
           CALL compute_entropy( entropy, f(1), nspin )
           entropy = entropy * nbsp
         END IF

         CALL vofrho( nfi, rhor(1,1), rhog(1,1), rhos(1,1), rhoc(1), tfirst, tlast, &
        &  ei1(1,1), ei2(1,1), ei3(1,1), irb(1,1), eigrb(1,1), sfac(1,1),           & 
        &  tau0(1,1), fion(1,1) )

         IF( tefield ) THEN
           CALL berry_energy( enb, enbi, bec, cm(:,:), fion ) 
           etot = etot + enb + enbi
         END IF
         IF( tefield2 ) THEN
           CALL berry_energy2( enb, enbi, bec, cm(:,:), fion )
           etot = etot + enb + enbi
         END IF

         CALL compute_stress( stress, detot, h, omega )

         if(iprsta.gt.2) CALL print_atomic_var( fion, na, nsp, ' fion ' )

         CALL newd( rhor, irb, eigrb, becsum, fion )
         !
         IF( force_pairing ) THEN
            !
            CALL runcp_uspp_force_pairing( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, cm, &
        &                 c0, ei_unp, fromscra = .TRUE. )
            ! lambda(nupdwn(1), nupdwn(1), 1) = ei_unp
            lambda(nupdwn(1), nupdwn(1), 2) = 0.d0 
            !
         ELSE
            !
#if defined __BGL
            CALL runcp_uspp_bgl( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, cm, &
        &                 c0(:,:), fromscra = .TRUE. )
#else
            CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, cm, &
        &                 c0(:,:), fromscra = .TRUE. )
#endif
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
            CALL ortho( eigr, c0, phi, ngw, lambda, SIZE(lambda,1), &
                        bigr, iter, ccc, bephi, becp, nbsp, nspin, nupdwn, iupdwn )
         else
            CALL gram( vkb, bec, nkb, c0, ngw, nbsp )
         endif
         !
         !
         if ( tfor .or. tprnfor ) CALL nlfl( bec, becdr, lambda, fion )

         if ( iprsta >= 3 ) CALL print_lambda( lambda, nbsp, 9, ccc )

         if ( tpre ) CALL nlfh( bec, dbec, lambda )
         !
         IF ( tortho ) THEN
            DO iss = 1, nspin_wfc
               CALL updatc( ccc, nbsp, lambda(:,:,iss), SIZE(lambda,1), phi, SIZE(phi,1), &
                            bephi, SIZE(bephi,1), becp, bec, c0, nupdwn(iss), iupdwn(iss) )
            END DO
         END IF
         !
         IF( force_pairing ) THEN
            !
            c0 ( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) =  c0( :, 1:nupdwn(2))
            phi( :, iupdwn(2):(iupdwn(2)+nupdwn(2)-1) ) = phi( :, 1:nupdwn(2))
            lambda(1:nupdwn(2), 1:nupdwn(2), 2) = lambda(1:nupdwn(2), 1:nupdwn(2), 1 )
            !
         ENDIF
         !
         CALL calbec ( nvb+1, nsp, eigr, c0, bec )

         if ( tpre ) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, cm, dbec, .true. )

         if(iprsta.ge.3) CALL dotcsc( eigr, c0, ngw, nbsp )
         !

         xnhp0=0.
         xnhpm=0.
         vnhp =0.
         fionm=0.
         CALL ions_vel( vels, taus, tausm, na, nsp, delt )
         xnhh0(:,:)=0.
         xnhhm(:,:)=0.
         vnhh (:,:) =0.
         velh (:,:)=(h(:,:)-hold(:,:))/delt
         !
         CALL elec_fakekine( ekincm, ema0bg, emass, c0, cm, ngw, nbsp, 1, delt )

         xnhe0=0.
         xnhem=0.
         vnhe =0.

         lambdam = lambda
         !
       ELSE 
          !
          c0 = cm
          !
       END IF

    ELSE

       edft%enl = nlrh( cm, ttforce, atoms%for, bec, becdr, eigr )
       !
       CALL rhoofr( 0, cm, cdesc, f, rhor, ht )
       !
       CALL vofrhos( ttprint, ttforce, tstress, rhor, atoms, &
                  vpot, bec, cm, cdesc, f, eigr, ei1, ei2, ei3, &
                  sfac, ht, edft )
       !
       IF( iprsta > 1 ) CALL debug_energies( edft )
       !
       IF ( tcarpar ) THEN
          !
          IF ( .NOT. force_pairing ) THEN
             !
             CALL runcp_ncpp( cm, cm, c0, vpot, vkb, f, bec, fccc, gam, fromscra = .TRUE. )
             !
          ELSE
             !
             c0 = cm
             !
          END IF
          !
          IF ( tortho .AND. ( .NOT. force_pairing ) ) THEN
             !
             CALL ortho( cm, c0, nupdwn, iupdwn, nspin )
             !
          ELSE
             !
             DO iss = 1, nspin_wfc
               !
               CALL gram( vkb, bec, nkb, c0(1,iupdwn(iss)), SIZE(c0,1), cdesc%nbt( iss ) )
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

       adum = 0.D0
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
