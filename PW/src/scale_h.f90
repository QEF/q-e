!
! Copyright (C) 2001-2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE scale_h
  !-----------------------------------------------------------------------
  !! When variable cell calculation are performed this routine scales the
  !! quantities needed in the calculation of the hamiltonian using the
  !! new and old cell parameters.
  !
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout
  USE cell_base,      ONLY : bg, omega, set_h_ainv, tpiba
  USE cellmd,         ONLY : at_old, omega_old
  USE constants,      ONLY : eps8
  USE gvect,          ONLY : g, gg, ngm
  USE klist,          ONLY : nks, xk, wk, ngk, igk_k, nkstot, qnorm
  USE control_flags,  ONLY : iverbosity
  USE start_k,        ONLY : nks_start, xk_start, nk1,nk2,nk3
  USE exx_base,       ONLY : exx_grid_init, exx_mp_init
  USE exx,            ONLY : exx_gvec_reinit
  USE xc_lib,         ONLY : xclib_dft_is
  USE rism_module,    ONLY : lrism, rism_reinit3d
  USE mp,             ONLY : mp_max
  USE mp_bands,       ONLY : intra_bgrp_comm
  USE mp_pools,       ONLY: inter_pool_comm
  USE atwfc_mod,      ONLY : scale_tab_atwfc, init_tab_atwfc
  USE beta_mod,       ONLY : scale_tab_beta, init_tab_beta
  USE qrad_mod,       ONLY : scale_tab_qrad, init_tab_qrad
  USE vloc_mod,       ONLY : scale_tab_vloc
  USE rhoc_mod,       ONLY : scale_tab_rhc
  USE rhoat_mod,      ONLY : scale_tab_rhoat
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, ik, ipol, ierr
  ! counters
  REAL(DP) :: gg_max, gk_max, gmax, kmax, gk2
  !
  ! scale the k points
  !
  CALL cryst_to_cart( nkstot, xk, at_old, - 1 )
  CALL cryst_to_cart( nkstot, xk, bg, + 1 )
  IF (nks_start>0) THEN
    CALL cryst_to_cart( nks_start, xk_start, at_old, - 1 )
    CALL cryst_to_cart( nks_start, xk_start, bg, + 1 )
  ENDIF
  !
  ! Print new k-points only if given in input and if not Gamma
  !
  IF ( nk1==0 .AND. nk2==0 .AND. nk3 == 0 .AND. &
       ( nkstot > 1 .OR. ABS(xk(1,1)**2+xk(2,1)**2+xk(3,1)**2) > eps8 ) ) THEN
     IF ( iverbosity > 0 .OR. nkstot < 100 ) THEN
        WRITE( stdout, '(5x,a)' ) 'NEW k-points:'
        DO ik = 1, nkstot
           WRITE( stdout, '(3f12.7,f12.7)') (xk(ipol,ik) , ipol=1,3), wk(ik)
        ENDDO
     ELSE
        WRITE( stdout, '(5x,a)' ) "NEW k-points: use verbosity='high' to print them"
     ENDIF
  ENDIF
  !
  ! scale the g vectors (as well as gg and gl arrays)
  !
  CALL cryst_to_cart( ngm, g, at_old, - 1 )
  CALL cryst_to_cart( ngm, g, bg, + 1 )
  gg_max = 0.0_dp
  !
  DO ig = 1, ngm
     gg (ig) = g(1,ig) * g(1,ig) + g(2,ig) * g(2,ig) + g(3,ig) * g(3,ig)
     gg_max = MAX(gg(ig), gg_max)
  ENDDO
  !$acc update device(g,gg)
  !
  CALL mp_max( gg_max, intra_bgrp_comm )
  gmax = SQRT(gg_max)*tpiba
  !
  ! gmax is the largest |G| actually needed in interpolation tables
  ! (gmax^2=ecutrho at the first step or at fixed cell)
  !
  IF ( xclib_dft_is('hybrid') ) THEN
     CALL exx_grid_init( reinit=.TRUE. )
     ! not sure next line is needed
     CALL exx_mp_init( )
     CALL exx_gvec_reinit( at_old )
     gmax = gmax + qnorm
     ! For USPP + hybrid, gmax is the largest |q+G| needed
     ! Note that qnorm is recomputed in exx_grid_init
  END IF
  !
  ! Now we need the largest |k+G| for beta and atwfc interpolation
  !
  gk_max = 0.0_dp
  DO ik = 1, nks
     DO ig = 1, ngk(ik)
        gk2 = (xk(1,ik) + g(1,igk_k(ig,ik)))**2 + &
              (xk(2,ik) + g(2,igk_k(ig,ik)))**2 + &
              (xk(3,ik) + g(3,igk_k(ig,ik)))**2
        gk_max = MAX(gk2, gk_max)
     END DO
  ENDDO
  CALL mp_max( gk_max, intra_bgrp_comm )
  CALL mp_max( gk_max, inter_pool_comm )
  kmax = SQRT(gk_max)*tpiba
  ! (kmax^2=ecutwfc at the first step or at fixed cell)
  ! for hybrid functionals, we need max |k+q+G|
  IF ( xclib_dft_is('hybrid') )  kmax = kmax + qnorm
  !
  ! Scale the interpolation tables with correct volume
  !
  CALL scale_tab_atwfc( omega_old/omega )
  CALL scale_tab_beta ( omega_old/omega )
  CALL scale_tab_rhc  ( omega_old/omega )
  CALL scale_tab_rhoat( omega_old/omega )
  CALL scale_tab_qrad ( omega_old/omega )
  CALL scale_tab_vloc ( omega_old/omega )
  !
  ! Check that interpolation tables are of sufficient size,
  ! re-allocate and re-compute if needed
  ! For tab_vloc, tab_rhc, tab_rhoc, this is done when used
  !
  WRITE( stdout, '(5x,"New effective cutoffs (rho, wfc):",2f8.2)' ) gmax**2,kmax**2
  CALL init_tab_qrad ( gmax, omega, intra_bgrp_comm, ierr)
  IF ( ierr == -1) WRITE( stdout, '(5x,"Interpolation table for Q(G) re-allocated")' ) 
  CALL init_tab_beta ( kmax, omega, intra_bgrp_comm, ierr)
  IF ( ierr == -1) WRITE( stdout, '(5x,"Interpolation table for beta(G) re-allocated")' ) 
  CALL init_tab_atwfc( kmax, omega, intra_bgrp_comm, ierr)
  IF ( ierr == -1) WRITE( stdout, '(5x,"Interpolation table for atomic wavefunctions re-allocated")' ) 
  !
  ! recalculate the local part of the pseudopotential
  !
  CALL init_vloc( )
  !
  ! for ts-vdw
  !
  CALL set_h_ainv()
  !
  ! for 3D-RISM
  !
  IF ( lrism ) CALL rism_reinit3d()
  !
  RETURN
  !
END SUBROUTINE scale_h

