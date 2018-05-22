!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine scale_h
  !-----------------------------------------------------------------------
  ! When variable cell calculation are performed this routine scales the
  ! quantities needed in the calculation of the hamiltonian using the
  ! new and old cell parameters.
  !
  USE kinds,      ONLY : dp
  USE io_global,  ONLY : stdout
  USE cell_base,  ONLY : bg, omega
  USE cellmd,     ONLY : at_old, omega_old
  USE constants,  ONLY : eps8
  USE gvect,      ONLY : g, gg, ngm
  USE klist,      ONLY : xk, wk, nkstot
  USE us,         ONLY : nqxq, qrad, tab, tab_at, dq
  USE control_flags, ONLY : iverbosity
  USE start_k,    ONLY : nks_start, xk_start, nk1,nk2,nk3
  USE exx_base,   ONLY : exx_grid_init, exx_mp_init
  USE exx,        ONLY : exx_gvec_reinit
  USE funct,      ONLY : dft_is_hybrid
  USE mp,         ONLY : mp_max
  USE mp_bands,   ONLY : intra_bgrp_comm
  !
  USE gvect_gpum, ONLY : using_g, using_g_d, using_gg, using_gg_d
  !
  implicit none
  !
  integer :: ig
  ! counter on G vectors

  integer  :: ik, ipol
  real(dp) :: gg_max
  !
  ! scale the k points
  !
  call cryst_to_cart (nkstot, xk, at_old, - 1)
  call cryst_to_cart (nkstot, xk, bg, + 1)
  IF (nks_start>0) THEN
    call cryst_to_cart (nks_start, xk_start, at_old, - 1)
    call cryst_to_cart (nks_start, xk_start, bg, + 1)
  ENDIF
  !
  ! Print new k-points only if given in input and if not Gamma
  !
  IF ( nk1==0 .AND. nk2==0 .AND. nk3 == 0 .AND. &
       ( nkstot > 1 .OR. ABS(xk(1,1)**2+xk(2,1)**2+xk(3,1)**2) > eps8 ) ) THEN
     IF ( iverbosity > 0 .OR. nkstot < 100 ) THEN
        WRITE( stdout, '(5x,a)' ) 'NEW k-points:'
        DO ik = 1, nkstot
           WRITE( stdout, '(3f12.7,f12.7)') (xk (ipol, ik) , ipol=1,3), wk (ik)
        ENDDO
     ELSE
        WRITE( stdout, '(5x,a)' ) "NEW k-points: use verbosity='high' to print them"
     ENDIF
  ENDIF
  !
  ! scale the g vectors (as well as gg and gl arrays)
  !
  call cryst_to_cart (ngm, g, at_old, - 1)
  call cryst_to_cart (ngm, g, bg, + 1)
  gg_max = 0.0_dp
  do ig = 1, ngm
     gg (ig) = g(1, ig) * g(1, ig) + g(2, ig) * g(2, ig) + g(3, ig) * g(3, ig)
     gg_max = max(gg(ig), gg_max)
  enddo

  CALL using_g(1); CALL using_gg(1)       ! g and gg are used almost only after
  CALL using_g_d(0); CALL using_gg_d(0) ! a single initialization.
  !                                                   This is a trick to avoid checking for sync everywhere.
  CALL mp_max (gg_max, intra_bgrp_comm)

  if(nqxq < int(sqrt(gg_max)/dq)+4) then
     call errore('scale_h', 'Not enough space allocated for radial FFT: '//&
                          'try restarting with a larger cell_factor.',1)
  endif
  !
  ! scale the non-local pseudopotential tables
  !
  tab(:,:,:) = tab(:,:,:) * sqrt (omega_old/omega)
  qrad(:,:,:,:) = qrad(:,:,:,:) * omega_old/omega
  tab_at(:,:,:) = tab_at(:,:,:) * sqrt (omega_old/omega)
  !
  ! recalculate the local part of the pseudopotential
  !
  call init_vloc ( )
  !
  ! for hybrid functionals
  !
  IF ( dft_is_hybrid() ) THEN
     CALL exx_grid_init( reinit=.true. )
     ! not sure next line is needed
     CALL exx_mp_init( )
     CALL exx_gvec_reinit( at_old )
  END IF
  !
  return
end subroutine scale_h

