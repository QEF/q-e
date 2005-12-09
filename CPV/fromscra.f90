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
  INTERFACE from_scratch
     MODULE PROCEDURE from_scratch_fpmd, from_scratch_cp
  END INTERFACE
  !
  CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE from_scratch_fpmd( rhoe, desc, cm, c0, cp, ce, cdesc, edesc, &
                                eigr, ei1, ei2, ei3, sfac, fi, ht, atoms, &
                                bec, becdr, vpot, edft )
    !------------------------------------------------------------------------
    !
    USE kinds,            ONLY : DP
    USE wave_types,       ONLY : wave_descriptor
    USE wave_functions,   ONLY : fixwave, wave_rand_init
    USE wave_base,        ONLY : wave_steepest
    USE charge_density,   ONLY : rhoofr
    USE cell_module,      only : boxdimensions
    USE cell_base,        ONLY : s_to_r
    USE electrons_module, ONLY : nspin, pmss, occn_init, occn_info
    USE ions_base,        ONLY : taui, cdmi, randpos
    USE mp,               ONLY : mp_end
    USE nl,               ONLY : nlrh_m
    USE energies,         ONLY : dft_energy_type, debug_energies
    USE potentials,       ONLY : vofrhos
    USE forces,           ONLY : dforce_all
    USE orthogonalize,    ONLY : ortho
    USE control_flags,    ONLY : tcarpar, tfor, thdyn, tortho, tpre, tranp, &
                                 force_pairing, iprsta, tprnfor, amprp
    USE charge_types,     ONLY : charge_descriptor
    USE time_step,        ONLY : delt
    USE runcp_module,     ONLY : runcp_ncpp
    use grid_dimensions,  only : nr1, nr2, nr3
    USE gvecp,            ONLY : ngm
    USE io_global,        ONLY : ionode, stdout
    USE parameters,       ONLY : nacx
    USE uspp,             ONLY : vkb, nkb
    !
    USE atoms_type_module,    ONLY : atoms_type
    USE phase_factors_module, ONLY : strucf, phfacs
    USE cp_electronic_mass,   ONLY : emass
    USE print_out_module,     ONLY : printout
    USE reciprocal_vectors,   ONLY : mill_l, gx
    !
    IMPLICIT NONE
    !
    COMPLEX(DP),             INTENT(OUT) :: eigr(:,:)
    COMPLEX(DP),             INTENT(OUT) :: ei1(:,:)
    COMPLEX(DP),             INTENT(OUT) :: ei2(:,:)
    COMPLEX(DP),             INTENT(OUT) :: ei3(:,:)
    COMPLEX(DP),             INTENT(OUT) :: sfac(:,:)
    REAL(DP),                INTENT(OUT) :: rhoe(:,:,:,:)
    REAL(DP),                INTENT(OUT) :: bec(:,:)
    REAL(DP),                INTENT(OUT) :: becdr(:,:,:)
    REAL(DP),                INTENT(OUT) :: fi(:,:,:)
    REAL(DP),                INTENT(OUT) :: vpot(:,:,:,:)
    TYPE(atoms_type)    ,    INTENT(OUT) :: atoms
    TYPE(dft_energy_type) ,  INTENT(OUT) :: edft
    TYPE(boxdimensions)  ,   INTENT(INOUT) :: ht 
    TYPE(charge_descriptor), INTENT(IN)    :: desc
    TYPE(wave_descriptor),   INTENT(IN)    :: cdesc, edesc
    COMPLEX(DP),             INTENT(INOUT) :: cm(:,:,:,:), c0(:,:,:,:)
    COMPLEX(DP),             INTENT(INOUT) :: cp(:,:,:,:), ce(:,:,:,:)
    !
    ! ... declare other variables
    !
    LOGICAL, PARAMETER :: ttprint = .TRUE.
    INTEGER, PARAMETER :: nfi = 0
    LOGICAL            :: ttforce
    LOGICAL            :: tstress
    COMPLEX(DP)        :: cgam(1,1,1)
    REAL(DP)           :: gam(1,1,1)
    INTEGER            :: ierr
    REAL(DP)           :: timepre, vdum
    REAL(DP)           :: s4, s5, cclock
    REAL(DP)           :: adum( nacx )
    REAL(DP)           :: hinv( 3, 3 )
    INTEGER            :: nspin_wfc
    INTEGER            :: iss
    !
    !
    ttforce = tfor  .or. tprnfor
    tstress = thdyn .or. tpre
    !
    IF ( ANY( tranp ) ) THEN
       !
       hinv = TRANSPOSE( ht%m1 )
       !
       CALL randpos( atoms%taus, atoms%na, &
                     atoms%nsp, tranp, amprp, hinv, atoms%mobile )
       !
       CALL s_to_r( atoms%taus, atoms%taur, atoms%na, atoms%nsp, ht%hmat )
       !
    END IF
    !
    CALL phfacs( ei1, ei2, ei3, eigr, mill_l, &
                 atoms%taus, nr1, nr2, nr3, atoms%nat )
    !
    CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngm )
    !
    ! ... initialize wave functions
    !
    nspin_wfc = nspin
    !
    IF( force_pairing ) nspin_wfc = 1
    !
    cm = 0.D0
    ce = 0.D0
    cp = 0.D0
    !
    DO iss = 1, nspin_wfc
       !
       CALL wave_rand_init( cm( :, :, 1, iss ) )
       !
    END DO
    !
    IF ( ionode ) &
       WRITE( stdout, fmt = '(//,3X, "Wave Initialization: random initial wave-functions" )' )
    !
    DO iss = 1, nspin_wfc
       !
       CALL gram( vkb, bec, nkb, cm(1,1,1,iss), SIZE(cm,1), cdesc%nbt( iss ) )
       !
    END DO
    !
    c0 = cm
    !
    ! ... initialize bands
    !
    CALL occn_init( fi )
    CALL occn_info( fi )
    !
    atoms%for = 0.D0
    !
    ! ... compute local form factors
    !
    CALL formf( .true. , edft%eself )
    !
    edft%enl = nlrh_m( cm, cdesc, ttforce, atoms, fi, bec, becdr, eigr )
    !
    CALL rhoofr( 0, cm, cdesc, fi, rhoe, desc, ht )
    !
    CALL vofrhos( ttprint, ttforce, tstress, rhoe, desc, atoms, &
                  vpot, bec, cm, cdesc, fi, eigr, ei1, ei2, ei3, &
                  sfac, timepre, ht, edft )
    !
    IF( iprsta > 1 ) CALL debug_energies( edft )
    !
    IF ( tcarpar ) THEN
       ! 
       IF ( .NOT. force_pairing ) THEN
          !
          CALL runcp_ncpp( cm, cm, c0, cdesc, vpot, eigr, fi, &
                           bec, vdum, gam, cgam, fromscra = .TRUE. )
          !
       ELSE
          !
          c0 = cm
          !
       END IF
       !
       IF ( tortho .AND. ( .NOT. force_pairing ) ) THEN
          !
          CALL ortho( cm, c0, cdesc, pmss, emass )
          !
       ELSE
          !
          DO iss = 1, nspin_wfc
            !
            CALL gram( vkb, bec, nkb, c0(1,1,1,iss), SIZE(c0,1), cdesc%nbt( iss ) )
            !
          END DO
          !
       END IF
       !
    ELSE
       !
       c0 = cm
       !
    END IF
    !
    adum = 0.D0
    !
    CALL printout( nfi, atoms, 0.0d0, 0.0d0, ttprint, ht, adum, adum, edft )
    !
    RETURN
    !
  END SUBROUTINE from_scratch_fpmd
  !
  !--------------------------------------------------------------------------
  SUBROUTINE from_scratch_cp( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst,   &
                              eself, fion, taub, irb, eigrb, b1, b2, b3, nfi,  &
                              rhog, rhor, rhos, rhoc, enl, ekin, stress, detot,&
                              enthal, etot, lambda, lambdam, lambdap, ema0bg,  &
                              dbec, delt, bephi, becp, velh, dt2bye, iforce,   &
                              fionm, xnhe0, xnhem, vnhe, ekincm )
    !--------------------------------------------------------------------------
    !
    USE kinds,                ONLY: DP
    USE control_flags,        ONLY : tranp, trane, trhor, iprsta, tpre,   &
                                     tzeroc, tzerop, tzeroe, tfor, thdyn, &
                                     lwf, tprnfor, tortho, amprp, ampre,  &
                                     tsde, ortho_eps, ortho_max
    USE ions_positions,       ONLY : taus, tau0, tausm, vels, velsm
    USE ions_base,            ONLY : na, nsp, randpos, zv, ions_vel, pmass
    USE cell_base,            ONLY : ainv, h, s_to_r, ibrav, omega, press, &
                                     hold, r_to_s, deth, wmass, iforceh,   &
                                     cell_force
    use electrons_base,       ONLY : nbsp
    USE energies,             ONLY : entropy
    USE uspp,                 ONLY : vkb, becsum, deeq, nkb
    USE wavefunctions_module, ONLY : c0, cm, phi => cp
    USE io_global,            ONLY : stdout
    USE cpr_subroutines,      ONLY : compute_stress, print_atomic_var, &
                                     print_lambda, elec_fakekine
    USE core,                 ONLY : nlcc_any
    USE gvecw,                ONLY : ngw
    USE reciprocal_vectors,   ONLY : gstart, mill_l
    USE gvecs,                ONLY : ngs
    USE wave_base,            ONLY : wave_steepest
    USE cvan,                 ONLY : nvb
    USE ions_nose,            ONLY : xnhp0,  xnhpm,  vnhp
    USE cell_nose,            ONLY : xnhh0, xnhhm,  vnhh
    USE cp_electronic_mass,   ONLY : emass
    USE efield_module,        ONLY : tefield, efield_berry_setup, &
                                     berry_energy, dforce_efield
    USE cg_module,            ONLY : tcg
    USE ensemble_dft,         ONLY : tens, compute_entropy
    USE runcp_module,         ONLY : runcp_uspp
    USE electrons_base,       ONLY : f, nspin
    USE phase_factors_module, ONLY : strucf
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
    REAL(DP)    :: lambda(:,:), lambdam(:,:), lambdap(:,:)
    REAL(DP)    :: ema0bg(:)
    REAL(DP)    :: dbec(:,:,:,:)
    REAL(DP)    :: delt
    REAL(DP)    :: bephi(:,:), becp(:,:)
    REAL(DP)    :: velh(:,:)
    REAL(DP)    :: dt2bye, xnhe0, xnhem, vnhe, ekincm
    !
    REAL(DP),    ALLOCATABLE :: emadt2(:), emaver(:)
    COMPLEX(DP), ALLOCATABLE :: c2(:), c3(:)
    REAL(DP)                 :: verl1, verl2
    REAL(DP)                 :: bigr
    INTEGER                  :: i, j, iter
    LOGICAL                  :: tlast = .FALSE.
    REAL(DP)                 :: fcell(3,3), ccc, enb, enbi, fccc
    !
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
    CALL phfac( tau0, ei1, ei2, ei3, eigr )
    !     
    CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
    !     
    CALL initbox ( tau0, taub, irb )
    !     
    CALL phbox( taub, eigrb )
    !     
    if( trane ) then
       !       
       !     random initialization
       !     
       CALL randin( 1, nbsp, gstart, ngw, ampre, cm )
      
    else 
       !       
       !     gaussian initialization
       !     
       ! CALL gausin( eigr, cm )  ! DEBUG to be check
      
    end if
    !
    ! ... prefor calculates vkb (used by gram)
    !
    CALL prefor( eigr, vkb )
    !
    CALL gram( vkb, bec, nkb, cm, ngw, nbsp )

    if( iprsta .ge. 3 ) CALL dotcsc( eigr, cm )

    hold = h
    velh = 0.0d0
    fion = 0.0d0
    tausm = taus
    vels  = 0.0d0
    lambdam = lambda
    !
    CALL formf( tfirst, eself )

    IF( tefield ) THEN
      CALL efield_berry_setup( eigr, tau0 )
    END IF

    IF( .NOT. tcg ) THEN

      CALL calbec ( 1, nsp, eigr, cm, bec )
      if (tpre) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, cm, dbec, .true. )

      CALL initbox ( tau0, taub, irb )
      CALL phbox( taub, eigrb )
!
      CALL rhoofr ( nfi, cm, irb, eigrb, bec, becsum, rhor, rhog, rhos, enl, ekin )
!
!     put core charge (if present) in rhoc(r)
!
      if ( nlcc_any ) CALL set_cc( irb, eigrb, rhoc )

      IF( tens ) THEN
        CALL compute_entropy( entropy, f(1), nspin )
        entropy = entropy * nbsp
      END IF

      CALL vofrho( nfi, rhor(1,1), rhog(1,1), rhos(1,1), rhoc(1), tfirst, tlast,                 &
        &  ei1(1,1), ei2(1,1), ei3(1,1), irb(1,1), eigrb(1,1), sfac(1,1), & 
        &  tau0(1,1), fion(1,1) )

      IF( tefield ) THEN
        CALL berry_energy( enb, enbi, bec, cm(:,:,1,1), fion ) 
        etot = etot + enb + enbi
      END IF

      CALL compute_stress( stress, detot, h, omega )

      if(iprsta.gt.2) CALL print_atomic_var( fion, na, nsp, ' fion ' )

      CALL newd( rhor, irb, eigrb, becsum, fion )
      CALL prefor( eigr, vkb )

!
      fccc = 0.0d0
      !
      CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, cm(:,:,1,1), c0(:,:,1,1), fromscra = .TRUE. )

!
!     nlfq needs deeq bec
!
      if( tfor .or. tprnfor ) CALL nlfq( cm, eigr, bec, becdr, fion )
!
!     calphi calculates phi
!     the electron mass rises with g**2
!
      CALL calphi( cm, ema0bg, bec, vkb, phi )


      if( tortho ) then
         CALL ortho( eigr, c0, phi, lambda, bigr, iter, ccc, ortho_eps, ortho_max, delt, bephi, becp )
      else
         CALL gram( vkb, bec, nkb, c0, ngw, nbsp )
      endif
!
!
      if ( tfor .or. tprnfor ) CALL nlfl( bec, becdr, lambda, fion )

      if ( iprsta >= 3 ) CALL print_lambda( lambda, nbsp, 9, ccc )

      if ( tpre ) CALL nlfh( bec, dbec, lambda )
!
      if ( tortho ) CALL updatc( ccc, lambda, phi, bephi, becp, bec, c0 )
      CALL calbec ( nvb+1, nsp, eigr, c0, bec )
      if ( tpre ) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, cm, dbec, .true. )

      if(iprsta.ge.3) CALL dotcsc(eigr,c0)
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
!     ======================================================
!     kinetic energy of the electrons
!     ======================================================

      CALL elec_fakekine( ekincm, ema0bg, emass, c0, cm, ngw, nbsp, delt )

      xnhe0=0.
      xnhem=0.
      vnhe =0.

      lambdam(:,:)=lambda(:,:)

    ELSE
      !
      ! ... Cojugate Gradient
      !
      c0 = cm
      !
    END IF
    !
    RETURN
    !
  END SUBROUTINE from_scratch_cp
  !
!=----------------------------------------------------------------------------=!
END MODULE from_scratch_module
!=----------------------------------------------------------------------------=!
