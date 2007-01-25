!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE compute_u_kq(ik, q)
  !-----------------------------------------------------------------------
  !
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE cell_base,            ONLY : at, bg, omega, tpiba, tpiba2
  USE wavefunctions_module, ONLY : evc
  USE klist,                ONLY : nkstot, xk
  USE wvfct,                ONLY : nbnd, nbndx, npwx, npw, igk, g2kin, current_k, btype
  USE gipaw_module
  USE gvect, only: g,ngm,ecutwfc,ngl,nrxx, nr1, nr2, nr3, nrx1, nrx2, nrx3
  USE uspp,   ONLY : nkb, vkb, okvan
  USE pwcom, ONLY : et
  USE becmod,              ONLY : becp
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE g_psi_mod, ONLY : h_diag, s_diag
  USE scf,           ONLY : vrs, vltot, vr
  USE noncollin_module,     ONLY : noncolin, npol
  USE complex_diis_module

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik
  REAL(DP), INTENT(IN) :: q(3)

  real(dp), allocatable :: etq(:)
  real(dp) :: xk_plus_q(3)
  integer :: ig, ntry, iter, notconv, ipol
  real(dp) :: v_of_0
  real(dp) :: avg_iter, dav_iter, cg_iter, ethr
  integer :: diis_iter
  logical :: lrot

  integer, parameter :: max_cg_iter = 200
  integer, parameter :: isolve = 0
  integer, parameter :: david = 4
  complex(dp) :: ZDOTC

  CALL start_clock( 'u_kq' )
  v_of_0 = SUM( vltot(1:nrxx) ) / DBLE( nr1 * nr2 * nr3 )
  CALL reduce( 1, v_of_0 )

  call gk_sort(xk(1,ik),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)
  current_k = ik
  if (lsda) current_spin = isk(ik)

  xk_plus_q(:) = xk(:,ik) + q(:)
  call init_us_2(npw,igk,xk_plus_q,vkb)
  
  ! ... read in wavefunctions from the SCF calculation at k
  ! ... and use it as a starting guess
  CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
  call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
  !g2kin(1:npw) = ( ( xk(1,ik) + g(1,igk(1:npw)) )**2 + &
  !                 ( xk(2,ik) + g(2,igk(1:npw)) )**2 + &
  !                 ( xk(3,ik) + g(3,igk(1:npw)) )**2 ) * tpiba2
  !call h_psi(npwx, npw, nbnd, evc, evq)
  !print*, ZDOTC(npw, evc(1,1), 1, evq(1,1), 1)*13.605692d0
  !print*, ZDOTC(npw, evc(1,2), 1, evq(1,2), 1)*13.605692d0
  ! initial guess
  evq = evc
 
  ! ... sets the kinetic energy
  g2kin(1:npw) = ( ( xk_plus_q(1) + g(1,igk(1:npw)) )**2 + &
                   ( xk_plus_q(2) + g(2,igk(1:npw)) )**2 + &
                   ( xk_plus_q(3) + g(3,igk(1:npw)) )**2 ) * tpiba2

  ALLOCATE( h_diag( npwx,npol ) )
  ALLOCATE( s_diag( npwx,npol ) )
  ALLOCATE(btype( nbnd, nkstot ) )

  nbndx = nbnd
  IF ( isolve == 0 ) nbndx = david * nbnd
  btype(:,:) = 1
  iter = 1
  allocate(etq(1:nbndx))
  ethr = conv_threshold

  !====================================================================
  ! Diagonalize
  !====================================================================
  IF ( isolve == 1 ) THEN
    ! ... Conjugate-Gradient diagonalization
    ! ... h_diag is the precondition matrix
    !!!call errore('compute_u_kq','CG not working properly?',-1)
    h_diag = 1.D0
    FORALL( ig = 1 : npw )
      h_diag(ig,:) = 1.D0 + g2kin(ig) + SQRT( 1.D0 + ( g2kin(ig) - 1.D0 )**2 )
    END FORALL

    ntry = 0
    CG_loop : DO
      IF ( ntry > 0 ) THEN
        CALL cinitcgg( npwx, npw, nbnd, nbnd, evq, evq, etq )
        avg_iter = avg_iter + 1.D0
      END IF

      CALL ccgdiagg( npwx, npw, nbnd, evq, etq, btype(1,ik), &
                     h_diag, ethr, max_cg_iter, .true., &
                     notconv, cg_iter )

      avg_iter = avg_iter + cg_iter
      ntry = ntry + 1                
      print*, etq(:)*13.605692d0

      ! ... exit condition
      IF ( test_exit_cond() ) EXIT  CG_loop
    END DO CG_loop


  ELSE IF ( isolve == 2 ) THEN
    ! ... RMM-DIIS method
    call errore('compute_u_kq','DIIS not working properly?',-1)
#if 0
    h_diag(1:npw) = g2kin(1:npw) + v_of_0
    CALL usnldiag( h_diag, s_diag )

    ntry = 0
    RMMDIIS_loop: DO
      CALL cdiisg( npw, npwx, nbnd, evq, etq, &
                   btype(1,ik), notconv, diis_iter, iter )

      avg_iter = avg_iter + diis_iter
      ntry = ntry + 1                

      ! ... exit condition
      IF ( test_exit_cond() ) EXIT  RMMDIIS_loop
    END DO RMMDIIS_loop
#endif

  ELSE IF (isolve == 0) THEN
    ! ... Davidson diagonalization
    ! ... h_diag are the diagonal matrix elements of the
    ! ... hamiltonian used in g_psi to evaluate the correction 
    ! ... to the trial eigenvectors
    do ipol = 1, npol
      h_diag(1:npw,ipol) = g2kin(1:npw) + v_of_0
    enddo
    CALL usnldiag( h_diag, s_diag )

    ntry = 0
    david_loop: DO
      lrot = ( iter == 1 )
      CALL cegterg( npw, npwx, nbnd, nbndx, evq, ethr, &
                    okvan, etq, btype(1,ik), notconv, &
                    lrot, dav_iter )

      avg_iter = avg_iter + dav_iter
      ntry = ntry + 1                

      ! ... exit condition
      IF ( test_exit_cond() ) EXIT david_loop                
    END DO david_loop 
    !!!print '(4F8.4)', etq(1:nbnd)*13.605692d0
  ELSE

    call errore('compute_u_kq', 'diagonalization method?', -1)

  END IF
  
  IF ( notconv > MAX( 5, nbnd / 4 ) ) THEN
    CALL errore( 'compute_u_kq', 'too many bands are not converged', 1 )
  END IF

  call stop_clock('u_kq')

  DEALLOCATE( h_diag )
  DEALLOCATE( s_diag )
  DEALLOCATE(btype )


  CONTAINS

    FUNCTION test_exit_cond()
    IMPLICIT NONE
    LOGICAL :: test_exit_cond
    test_exit_cond =  ( ntry > 100 ) .OR.  ( notconv <= 0 )
    END FUNCTION test_exit_cond


END SUBROUTINE compute_u_kq
