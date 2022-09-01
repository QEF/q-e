!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
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
  USE gvect,          ONLY : g, gg, ngm, g_d, gg_d
  USE klist,          ONLY : xk, wk, nkstot
  USE uspp_data,      ONLY : nqxq, dq, scale_uspp_data
  USE control_flags,  ONLY : iverbosity
  USE start_k,        ONLY : nks_start, xk_start, nk1,nk2,nk3
  USE exx_base,       ONLY : exx_grid_init, exx_mp_init
  USE exx,            ONLY : exx_gvec_reinit
  USE xc_lib,         ONLY : xclib_dft_is
  USE rism_module,    ONLY : lrism, rism_reinit3d
  USE mp,             ONLY : mp_max
  USE mp_bands,       ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ig
  ! counter on G vectors
  INTEGER  :: ik, ipol
  REAL(DP) :: gg_max
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
#if defined(__CUDA)
  ! update GPU copies of variables as well
  g_d  = g
  gg_d = gg
#endif
  !$acc update device(g,gg)
  !
  CALL mp_max( gg_max, intra_bgrp_comm )
  !
  IF (nqxq < INT(SQRT(gg_max)*tpiba/dq)+4) THEN
     CALL errore( 'scale_h', 'Not enough space allocated for radial FFT: '//&
                             'try restarting with a larger cell_factor.', 1 )
  ENDIF
  !
  ! scale the non-local pseudopotential tables
  !
  call scale_uspp_data( omega_old/omega )

  !
  ! recalculate the local part of the pseudopotential
  !
  CALL init_vloc( )
  !
  ! for hybrid functionals
  !
  IF ( xclib_dft_is('hybrid') ) THEN
     CALL exx_grid_init( reinit=.TRUE. )
     ! not sure next line is needed
     CALL exx_mp_init( )
     CALL exx_gvec_reinit( at_old )
  ENDIF
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

