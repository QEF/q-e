!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"

!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
  MODULE from_scratch_module
!=----------------------------------------------------------------------------=!

        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: from_scratch

        INTERFACE from_scratch
          MODULE PROCEDURE from_scratch_fpmd, from_scratch_cp
        END INTERFACE

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE from_scratch_fpmd( kp, ps, rhoe, desc, cm, c0, cdesc, &
        eigr, ei1, ei2, ei3, sfac, fi, ht, atoms, fnl, vpot, edft )

!  (describe briefly what this routine does...)
!  ----------------------------------------------
!  END manual


! ... declare modules
      USE kinds
      USE cp_types, ONLY: pseudo
      USE wave_types, ONLY: wave_descriptor
      USE atoms_type_module, ONLY: atoms_type
      USE wave_functions, ONLY: gram, fixwave
      USE wave_base, ONLY: wave_steepest
      USE charge_density, ONLY: rhoofr
      USE phase_factors_module, ONLY: strucf, phfacs
      USE cell_module, only: boxdimensions
      USE electrons_module, ONLY: nspin, pmss
      USE cp_electronic_mass, ONLY: emass
      USE ions_base, ONLY: taui, cdmi
      USE ions_module, ONLY: set_reference_positions
      USE mp, ONLY: mp_end
      USE nl, ONLY: nlrh_m
      USE energies, ONLY: dft_energy_type, debug_energies
      USE potentials, ONLY: vofrhos
      USE forces, ONLY: dforce_all
      USE orthogonalize, ONLY: ortho
      USE brillouin, ONLY: kpoints
      USE pseudo_projector, ONLY: projector
      USE control_flags, ONLY: tcarpar, tfor, thdyn, tortho, force_pairing
      USE charge_types, ONLY: charge_descriptor
      USE time_step, ONLY: delt
      USE runcp_module, ONLY: runcp_ncpp
      use grid_dimensions,    only: nr1, nr2, nr3
      USE reciprocal_vectors, ONLY: mill_l
      USE gvecp, ONLY: ngm

      IMPLICIT NONE

! ... declare subroutine arguments

      TYPE (atoms_type) :: atoms
      COMPLEX(dbl) :: eigr(:,:)
      COMPLEX(dbl) :: ei1(:,:)
      COMPLEX(dbl) :: ei2(:,:)
      COMPLEX(dbl) :: ei3(:,:)
      TYPE (kpoints) :: kp 
      REAL(dbl) :: rhoe(:,:,:,:)
      COMPLEX(dbl) :: sfac(:,:)
      TYPE (charge_descriptor), INTENT(IN) :: desc
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      COMPLEX(dbl), INTENT(INOUT) :: cm(:,:,:,:), c0(:,:,:,:)
      TYPE (boxdimensions) :: ht 
      TYPE (pseudo) :: ps
      REAL(dbl) :: fi(:,:,:)
      TYPE (projector) :: fnl(:,:)
      REAL (dbl)    ::  vpot(:,:,:,:)
      TYPE (dft_energy_type) :: edft


! ... declare other variables

      LOGICAL, PARAMETER :: ttforce = .TRUE.
      LOGICAL, PARAMETER :: ttprint = .TRUE.
      INTEGER, PARAMETER :: nfi = 0

      COMPLEX(dbl) :: cgam(1,1,1)
      REAL (dbl)   :: gam(1,1,1)
      INTEGER :: ierr
      REAL(dbl)  timepre, vdum

!  end of declarations
!  ----------------------------------------------


      atoms%for = 0.0d0
      CALL phfacs( ei1, ei2, ei3, eigr, mill_l, atoms%taus, nr1, nr2, nr3, atoms%nat )
      CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngm )
      edft%enl = nlrh_m(cm, cdesc, ttforce, atoms, fi, kp, fnl, ps%wsg, ps%wnl, eigr)
      CALL rhoofr( 0, kp, cm, cdesc, fi, rhoe, desc, ht)
      CALL vofrhos(ttprint, rhoe, desc, tfor, thdyn, ttforce, atoms, &
           kp, fnl, vpot, ps, cm, cdesc, fi, eigr, ei1, ei2, ei3, sfac, timepre, ht, edft)

      CALL debug_energies( edft ) ! DEBUG

      IF( tcarpar ) THEN
        
        IF( .NOT. force_pairing ) THEN

          CALL runcp_ncpp( cm, cm, c0, cdesc, kp, ps, vpot, eigr, &
             fi, fnl, vdum, gam, cgam, fromscra = .TRUE. )

        ELSE

          c0 = cm

        END IF

        IF( tortho .AND. ( .NOT. force_pairing ) ) THEN
          CALL ortho( cm, c0, cdesc, pmss, emass )
        ELSE
          CALL gram( c0, cdesc )
        END IF

      ELSE

        c0 = cm

      END IF

      CALL set_reference_positions(cdmi, taui, atoms, ht)

   RETURN
   END SUBROUTINE from_scratch_fpmd

!=----------------------------------------------------------------------------=!

SUBROUTINE from_scratch_cp( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst, eself, fion, &
      taub, irb, eigrb, b1, b2, b3, nfi, rhog, rhor, rhos, rhoc, enl, ekin, stress,  &
      detot, enthal, etot, lambda, lambdam, lambdap, ema0bg, dbec, delt,  &
      bephi, becp, velh, dt2bye, iforce, fionm, nbeg, xnhe0, xnhem, vnhe, ekincm )

    USE control_flags, ONLY: tranp, trane, trhor, iprsta, tpre, tzeroc 
    USE control_flags, ONLY: tzerop, tzeroe, tfor, thdyn, lwf, tprnfor, tortho
    USE control_flags, ONLY: amprp, taurdr, ampre, tsde, ortho_eps, ortho_max
    USE ions_positions, ONLY: taus, tau0, tausm, vels, velsm
    USE ions_base, ONLY: na, nsp, randpos, zv, ions_vel, pmass
    USE cell_base, ONLY: ainv, h, s_to_r, ibrav, omega, press, hold, r_to_s, deth
    USE cell_base, ONLY: wmass, iforceh, cell_force
    use electrons_base, only: n => nbsp
    USE energies, ONLY: entropy
    USE uspp, ONLY: betae => vkb, rhovan => becsum, deeq
    USE wavefunctions_module, ONLY: c0, cm, phi => cp
    USE io_global, ONLY: stdout
    USE cpr_subroutines, ONLY: compute_stress, print_atomic_var, print_lambda, elec_fakekine
    USE core, ONLY: nlcc_any
    USE gvecw, ONLY: ngw
    USE reciprocal_vectors, ONLY: gstart, mill_l
    USE gvecs, ONLY: ngs
    USE wave_base, ONLY: wave_steepest
    USE cvan, ONLY: nvb
    USE ions_nose, ONLY: xnhp0,  xnhpm,  vnhp
    USE cell_nose, ONLY: xnhh0, xnhhm,  vnhh
    USE cp_electronic_mass, ONLY: emass, emaec => emass_cutoff
    USE efield_module, ONLY: tefield, efield_berry_setup, berry_energy, dforce_efield
    USE cg_module, ONLY: tcg
    USE ensemble_dft, ONLY: tens, compute_entropy
    USE runcp_module, ONLY: runcp_uspp
    USE electrons_base, ONLY: f, nspin
    USE phase_factors_module, ONLY: strucf

    COMPLEX(kind=8) :: eigr(:,:), ei1(:,:),  ei2(:,:),  ei3(:,:)
    COMPLEX(kind=8) :: eigrb(:,:)
    REAL(kind=8) :: bec(:,:), fion(:,:), becdr(:,:,:), fionm(:,:)
    REAL(kind=8) :: eself
    REAL(kind=8) :: taub(:,:)
    REAL(kind=8) :: b1(:), b2(:), b3(:)
    INTEGER :: irb(:,:)
    INTEGER :: nfi, iforce(:,:), nbeg
    LOGICAL :: tfirst
    COMPLEX(kind=8) :: sfac(:,:)
    COMPLEX(kind=8) :: rhog(:,:)
    REAL(kind=8) :: rhor(:,:), rhos(:,:), rhoc(:), enl, ekin
    REAL(kind=8) :: stress(:,:), detot(:,:), enthal, etot
    REAL(kind=8) :: lambda(:,:), lambdam(:,:), lambdap(:,:)
    REAL(kind=8) :: ema0bg(:)
    REAL(kind=8) :: dbec(:,:,:,:)
    REAL(kind=8) :: delt
    REAL(kind=8) :: bephi(:,:), becp(:,:)
    REAL(kind=8) :: velh(:,:)
    REAL(kind=8) :: dt2bye, xnhe0, xnhem, vnhe, ekincm


    REAL(kind=8), ALLOCATABLE :: emadt2(:), emaver(:)
    COMPLEX(kind=8), ALLOCATABLE :: c2(:), c3(:)
    REAL(kind=8) :: verl1, verl2
    REAL(kind=8) :: bigr
    INTEGER :: i, j, iter
    LOGICAL :: tlast = .FALSE.
    REAL(kind=8) :: fcell(3,3), ccc, dt2hbe, enb, enbi, fccc

    !
    !

    ! WRITE(6,*) 'DEBUG running from_scratch'

    dt2hbe = 0.5d0 * dt2bye

    ! Input positions read from input file and stored in tau0

    IF( taurdr ) THEN
      call r_to_s( tau0, taus, na, nsp, h )
    END IF

    IF( ANY( tranp( 1:nsp ) ) ) THEN
      call invmat( 3, h, ainv, deth )
      call randpos(taus, na, nsp, tranp, amprp, ainv, iforce )
      call s_to_r( taus, tau0, na, nsp, h )
    END IF

!
    call phfac( tau0, ei1, ei2, ei3, eigr )
!     
    call initbox ( tau0, taub, irb )
!     
    call phbox( taub, eigrb )
!     
    if( trane ) then
      !       
      !     random initialization
      !     
      call randin( 1, n, gstart, ngw, ampre, cm )
      
    else 
      !       
      !     gaussian initialization
      !     
      ! call gausin( eigr, cm )  ! DEBUG to be check
      
    end if

    call prefor( eigr, betae )   !     prefor calculates betae (used by gram)
    call gram( betae, bec, cm )

    if( iprsta .ge. 3 ) call dotcsc( eigr, cm )

    hold = h
    velh = 0.0d0
    fion = 0.0d0
    tausm = taus
    vels  = 0.0d0
    lambdam = lambda
!     
    call phfac( tau0, ei1, ei2, ei3, eigr )
    call strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
    call formf( tfirst, eself )

    IF( tefield ) THEN
      CALL efield_berry_setup( eigr, tau0 )
    END IF



    IF( .NOT. tcg ) THEN

      call calbec ( 1, nsp, eigr, cm, bec )
      if (tpre) call caldbec( 1, nsp, eigr, cm )

      call initbox ( tau0, taub, irb )
      call phbox( taub, eigrb )
!
      call rhoofr ( nfi, cm, irb, eigrb, bec, rhovan, rhor, rhog, rhos, enl, ekin )
!
!     put core charge (if present) in rhoc(r)
!
      if ( nlcc_any ) call set_cc( irb, eigrb, rhoc )

      IF( tens ) THEN
        CALL compute_entropy( entropy, f(1), nspin )
        entropy = entropy * n
      END IF

      call vofrho( nfi, rhor(1,1), rhog(1,1), rhos(1,1), rhoc(1), tfirst, tlast,                 &
        &  ei1(1,1), ei2(1,1), ei3(1,1), irb(1,1), eigrb(1,1), sfac(1,1), & 
        &  tau0(1,1), fion(1,1) )

      IF( tefield ) THEN
        CALL berry_energy( enb, enbi, bec, cm(:,:,1,1), fion ) 
        etot = etot + enb + enbi
      END IF

      call compute_stress( stress, detot, h, omega )

      if(iprsta.gt.2) call print_atomic_var( fion, na, nsp, ' fion ' )

      call newd( rhor, irb, eigrb, rhovan, fion )
      call prefor( eigr, betae )

!
      fccc = 0.0d0
      !
      CALL runcp_uspp( nfi, fccc, ccc, ema0bg, dt2bye, rhos, bec, cm(:,:,1,1), c0(:,:,1,1), fromscra = .TRUE. )

!
!     nlfq needs deeq bec
!
      if( tfor .or. tprnfor ) call nlfq( cm, eigr, bec, becdr, fion )
!
!     calphi calculates phi
!     the electron mass rises with g**2
!
      call calphi( cm, ema0bg, bec, betae, phi )


      if( tortho ) then
         call ortho( eigr, c0, phi, lambda, bigr, iter, ccc, ortho_eps, ortho_max, delt, bephi, becp )
      else
         call gram( betae, bec, c0 )
      endif
!
!
      if ( tfor .or. tprnfor ) call nlfl( bec, becdr, lambda, fion )

      if ( iprsta >= 3 ) call print_lambda( lambda, n, 9, ccc )

      if ( tpre ) call nlfh( bec, dbec, lambda )
!
      if ( tortho ) call updatc( ccc, lambda, phi, bephi, becp, bec, c0 )
      call calbec ( nvb+1, nsp, eigr, c0, bec )
      if ( tpre ) call caldbec( 1, nsp, eigr, cm )

      if(iprsta.ge.3) call dotcsc(eigr,c0)
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

      call elec_fakekine( ekincm, ema0bg, emass, c0, cm, ngw, n, delt )

      xnhe0=0.
      xnhem=0.
      vnhe =0.

      lambdam(:,:)=lambda(:,:)


    else

      !
      !  Cojugate Gradient

      c0 = cm

    end if

    return
  END SUBROUTINE from_scratch_cp


!=----------------------------------------------------------------------------=!
      END MODULE from_scratch_module
!=----------------------------------------------------------------------------=!
