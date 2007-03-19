
! Copyright (C) 2001-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE compute_u_kq(ik, q)
  !----------------------------------------------------------------------------
  !
  ! ... diagonalize at k+q
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunigk, nwordatwfc, iunsat, iunwfc, &
                                   nwordwfc
  USE klist,                ONLY : nkstot, nks, xk, ngk
  USE uspp,                 ONLY : vkb, nkb
  USE gvect,                ONLY : g, nrxx, nr1, nr2, nr3  
  USE wvfct,                ONLY : et, nbnd, npwx, igk, npw, g2kin, &
                                   current_k, nbndx, btype
  USE control_flags,        ONLY : ethr, isolve, io_level, lscf
  USE ldaU,                 ONLY : lda_plus_u, swfcatom
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions_module, ONLY : evc  
  USE gvect,                ONLY : g, ngm, ecutwfc, ngl, nrxx, &
                                   nr1, nr2, nr3, nrx1, nrx2, nrx3
  USE cell_base,            ONLY : at, bg, omega, tpiba, tpiba2
  USE bp,                   ONLY : lelfield
  USE control_flags,        ONLY : iverbosity
  USE becmod,               ONLY : becp
  USE control_flags,        ONLY : ethr, isolve, max_cg_iter
  USE buffers
  USE gipaw_module
  IMPLICIT NONE
  INTEGER :: ik, iter       ! k-point, current iterations
  REAL(DP) :: q(3)          ! q-vector
  REAL(DP) :: avg_iter
  INTEGER :: ig
  REAL(DP) :: xkold(3)
  REAL(DP), allocatable :: et_old(:,:)

  CALL start_clock( 'c_bands' )

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  iter = 1
  max_cg_iter = 200
  isolve = 1  ! CG
  nbndx = nbnd
  ethr = conv_threshold
  if (allocated(btype)) deallocate(btype)
  allocate(btype(nbnd,nkstot))
  btype(:,:) = 1
  lscf = .false.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! save eigenvalues
  allocate( et_old(nbnd,nkstot) )
  et_old = et

  !! debug 
  !!WRITE(stdout, '(5X,"compute_u_kq: q = (",F10.4,",",F10.4,",",F10.4,")")') q
  !!WRITE(stdout, '(5X,"  isolve = ", I2)') isolve

  avg_iter = 0.D0

  current_k = ik
  IF ( lsda ) current_spin = isk(ik)

  ! same sorting of G-vector at k+q
  call gk_sort(xk(1,ik),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)

  ! set the k-point
  xkold(:) = xk(:,ik)
  xk(:,ik) = xk(:,ik) + q(:)
  g2kin(1:npw) = ( ( xk(1,ik) + g(1,igk(1:npw)) )**2 + &
                   ( xk(2,ik) + g(2,igk(1:npw)) )**2 + &
                   ( xk(3,ik) + g(3,igk(1:npw)) )**2 ) * tpiba2

  ! various initializations
  IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )

  ! read in wavefunctions from the previous iteration
  !!IF ( nks > 1 .OR. (io_level > 1) .OR. lelfield ) &
  CALL get_buffer( evc, nwordwfc, iunwfc, ik)

  ! Needed for LDA+U
  IF ( lda_plus_u ) CALL davcio( swfcatom, nwordatwfc, iunsat, ik, -1 )

  ! diagonalization of bands for k-point ik
  call diag_bands ( iter, ik, avg_iter )

  call poolreduce( 1, avg_iter )
  avg_iter = avg_iter / nkstot

  !! debug
  !!write( stdout, &
  !!     '( 5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1 )' ) &
  !!     ethr, avg_iter

  ! restore the k-point and eigenvalues
  !print*, et(1:nbnd,ik)*13.6056982d0
  xk(:,ik) = xkold(:)
  et = et_old
  deallocate(et_old)

  ! restore wavefunctions
  evq = evc
  CALL get_buffer( evc, nwordwfc, iunwfc, ik)

  CALL stop_clock( 'c_bands' )  
  RETURN

END SUBROUTINE compute_u_kq
