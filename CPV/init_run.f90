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
SUBROUTINE init_run()
  !----------------------------------------------------------------------------
  !
  ! ... this routine initialise the cp code and allocates (calling the
  ! ... appropriate routines) the memory
  !
  USE kinds,                    ONLY : dbl
  USE control_flags,            ONLY : nfi, nbeg, nomore, lwf, iprsta, &
                                       ndr, tfor, tprnfor, tpre
  USE cp_electronic_mass,       ONLY : emass, emass_cutoff
  USE ions_base,                ONLY : na, nax, nat, nsp, iforce, pmass, &
                                       fion, fionm, cdmi
  USE ions_positions,           ONLY : tau0, taum, taup, taus, tausm, tausp, &
                                       vels, velsm, velsp
  USE gvecw,                    ONLY : ngw, ecutw
  USE gvecb,                    ONLY : ngb
  USE gvecs,                    ONLY : ngs
  USE gvecp,                    ONLY : ngm
  USE grid_dimensions,          ONLY : nnrx, nr1, nr2, nr3
  USE electrons_base,           ONLY : nspin, nbsp, nbspx, nupdwn, f
  USE uspp,                     ONLY : nkb, vkb, deeq, becsum
  USE core,                     ONLY : nlcc_any, rhoc
  USE smooth_grid_dimensions,   ONLY : nnrsx
  USE wavefunctions_module,     ONLY : c0, cm, cp
  USE cdvan,                    ONLY : dbec, drhovan
  USE derho,                    ONLY : drhor, drhog
  USE ensemble_dft,             ONLY : tens, z0
  USE cg_module,                ONLY : tcg
  USE electrons_base,           ONLY : nudx
  USE parameters,               ONLY : natx
  USE efield_module,            ONLY : tefield
  USE uspp_param,               ONLY : nhm
  USE ions_nose,                ONLY : xnhp0, xnhpm, vnhp, nhpcl
  USE cell_base,                ONLY : h, hold, hnew, velh, tpiba2, ibrav, &
                                       alat, celldm, a1, a2, a3, b1, b2, b3
  USE cp_main_variables,        ONLY : lambda, lambdam, lambdap, ema0bg, bec,  &
                                       becdr, sfac, eigr, ei1, ei2, ei3, taub, &
                                       irb, eigrb, rhog, rhos, rhor, bephi,    &
                                       becp, acc
  USE gvecw,                    ONLY : ggp
  USE energies,                 ONLY : eself, enl, ekin, etot, enthal, ekincm
  USE stre,                     ONLY : stress
  USE dener,                    ONLY : detot
  USE time_step,                ONLY : dt2, delt, tps
  USE electrons_nose,           ONLY : xnhe0, xnhem, vnhe, fccc
  USE cell_nose,                ONLY : xnhh0, xnhhm, vnhh
  USE gvecp,                    ONLY : ecutp
  USE metagga,                  ONLY : ismeta, crosstaus, dkedtaus, gradwfc
  !
  USE cp_main_variables,        ONLY : allocate_mainvar
  USE efcalc,                   ONLY : clear_nbeg
  USE cpr_subroutines,          ONLY : print_atomic_var
  USE local_pseudo,             ONLY : allocate_local_pseudo
  USE cp_electronic_mass,       ONLY : emass_precond
  USE wannier_subroutines,      ONLY : wannier_startup
  USE from_scratch_module,      ONLY : from_scratch
  USE from_restart_module,      ONLY : from_restart  
  USE restart_file,             ONLY : readfile
  USE ions_base,                ONLY : ions_cofmass
  USE ensemble_dft,             ONLY : id_matrix_init, allocate_ensemble_dft
  USE efield_module,            ONLY : allocate_efield
  USE cg_module,                ONLY : allocate_cg
  USE wannier_module,           ONLY : allocate_wannier  
  !
  IMPLICIT NONE
  !
  !
  CALL start_clock( 'initialize' )
  !
  ! ... initialize g-vectors, fft grids
  !
  CALL init_dimensions()
  !
  CALL init_geometry()
  !
  IF ( lwf ) CALL clear_nbeg( nbeg )
  !
  ! ... more initialization requiring atomic positions
  !
  nax = MAXVAL( na(1:nsp) )
  !
  IF ( iprsta > 1 ) CALL print_atomic_var( tau0, na, nsp, ' tau0 ' )
  !
  !=======================================================================
  !     allocate and initialize nonlocal potentials
  !=======================================================================
  !
  CALL nlinit()
  !
  !=======================================================================
  !     allocation of all arrays not already allocated in init and nlinit
  !=======================================================================
  !
  CALL allocate_mainvar( ngw, ngb, ngs, ngm, nr1, nr2, nr3, nnrx, &
                         nnrsx, nat, nax, nsp, nspin, nbsp, nbspx, nkb )
  !
  IF ( nlcc_any ) THEN
     !
     ALLOCATE( rhoc( nnrx ) )
     !
  ELSE
     !
     ! ... dummy allocation required because this array appears in the
     ! ... list of arguments of some routines
     !
     ALLOCATE( rhoc( 1 ) )
     !
  END IF
  !
  ALLOCATE( c0( ngw, nbspx, 1, 1 ) )
  ALLOCATE( cm( ngw, nbspx, 1, 1 ) )
  ALLOCATE( cp( ngw, nbspx, 1, 1 ) )
  !
  CALL allocate_local_pseudo( ngs, nsp )
  !
  ALLOCATE( vkb( ngw, nkb ) )
  ALLOCATE( deeq( nhm, nhm, nat, nspin ) )
  ALLOCATE( dbec( nkb, nbsp, 3, 3 ) )
  ALLOCATE( drhog( ngm,  nspin, 3, 3 ) )
  ALLOCATE( drhor( nnrx, nspin, 3, 3 ) )
  ALLOCATE( becsum(  nhm*(nhm+1)/2, nat, nspin ) )
  ALLOCATE( drhovan( nhm*(nhm+1)/2, nat, nspin, 3, 3 ) )
  !
  IF ( ismeta .AND. tens ) &
     CALL errore( 'cprmain ', 'ensemble_dft not implimented for metaGGA', 1 )
  !
  IF ( ismeta .AND. tpre ) THEN
     !
     ALLOCATE( crosstaus( nnrsx, 6, nspin ) )
     ALLOCATE( dkedtaus(  nnrsx, 3, 3, nspin ) )
     ALLOCATE( gradwfc(   nnrsx, 3 ) )
     !
  END IF
  !
  IF ( lwf ) CALL allocate_wannier( nbsp, nnrsx, nspin, ngm )
  !
  IF ( tens .OR. tcg) &
     CALL allocate_ensemble_dft( nkb, nbsp, ngw, nudx, &
                                 nspin, nbspx, nnrsx, natx )
  !
  IF( tcg ) CALL allocate_cg( ngw, nbspx )
  !
  IF( tefield ) CALL allocate_efield( ngw, nbspx, nhm, nax, nsp )
  !
  deeq(:,:,:,:) = 0.D0
  !
  !=======================================================================
  !
  taum  = tau0
  taup  = 0.D0
  tausm = taus
  tausp = 0.D0
  vels  = 0.D0
  velsm = 0.D0
  velsp = 0.D0
  !
  hnew = h
  !
  lambda(:,:) = 0.D0
  cm(:,:,1,1) = ( 0.D0, 0.D0 )
  c0(:,:,1,1) = ( 0.D0, 0.D0 )
  !
  IF ( tens ) CALL id_matrix_init( nupdwn, nspin )
  !
  CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emass_cutoff )
  !
  IF ( lwf ) CALL wannier_startup( ibrav, alat, a1, a2, a3, b1, b2, b3 )
  !
  IF ( nbeg < -1 ) THEN
     !
     !======================================================================
     !     Initialize from scratch nbeg = -2 or nbeg = -3
     !======================================================================
     !
     nfi = 0
     !
     CALL from_scratch( sfac, eigr, ei1, ei2, ei3, bec, becdr, .TRUE.,    &
                        eself, fion, taub, irb, eigrb, b1, b2, b3, nfi,   &
                        rhog, rhor, rhos, rhoc, enl, ekin, stress, detot, &
                        enthal, etot, lambda, lambdam, lambdap, ema0bg,   &
                        dbec, delt, bephi, becp, velh, dt2/emass, iforce, &
                        fionm, nbeg, xnhe0, xnhem, vnhe, ekincm )
     !
  ELSE
     !
     !======================================================================
     !     nbeg = -1, nbeg = 0, nbeg = 1 or nbeg = 2
     !======================================================================
     !
     CALL readfile( 1, ndr, h, hold, nfi, c0(:,:,1,1), cm(:,:,1,1), taus,   &
                    tausm, vels, velsm, acc, lambda, lambdam, xnhe0, xnhem, &
                    vnhe, xnhp0, xnhpm, vnhp, nhpcl, ekincm, xnhh0, xnhhm,  &
                    vnhh, velh, ecutp, ecutw, delt, pmass, ibrav, celldm,   &
                    fion, tps, z0, f )
     !
     CALL from_restart( sfac, eigr, ei1, ei2, ei3, bec, becdr, .TRUE., fion, &
                        taub, irb, eigrb, b1, b2, b3, nfi, rhog, rhor, rhos, &
                        rhoc, stress, detot, enthal, lambda, lambdam,        &
                        lambdap, ema0bg, dbec, bephi, becp, velh, dt2/emass, &
                        fionm, ekincm )
     !
  END IF
  !
  !=======================================================================
  !     restart with new averages and nfi=0
  !=======================================================================
  !
  ! ... Fix. Center of Mass - M.S
  !
  IF ( lwf ) CALL ions_cofmass( tau0, pmass, na, nsp, cdmi )
  !
  IF ( nbeg <= 0 ) THEN
     !
     acc = 0.D0
     nfi = 0
     !
  END IF
  !
  IF ( .NOT. tfor .AND. .NOT. tprnfor ) fion(:,:) = 0.D0
  !
  IF ( .NOT. tpre ) stress = 0.D0
  !         
  fccc = 1.D0
  !
  nomore = nomore + nfi
  !
  CALL ions_cofmass( taus, pmass, na, nsp, cdmi )
  !
  CALL stop_clock( 'initialize' )
  !
  RETURN
  !
END SUBROUTINE init_run
