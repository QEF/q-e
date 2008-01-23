
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
  USE constants,            ONLY : RytoeV, tpi
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunigk, nwordatwfc, iunsat, iunwfc, &
                                   nwordwfc
  USE klist,                ONLY : nkstot, nks, xk, ngk
  USE uspp,                 ONLY : vkb, nkb
  USE wvfct,                ONLY : et, nbnd, npwx, igk, npw, g2kin, &
                                   current_k, nbndx, btype
  USE control_flags,        ONLY : ethr, io_level, lscf, istep, max_cg_iter
  USE control_flags,        ONLY : cntrl_isolve => isolve
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
  USE random_numbers,       ONLY : rndm
  USE buffers
  USE gipaw_module
  IMPLICIT NONE
  INTEGER :: ik, iter       ! k-point, current iterations
  REAL(DP) :: q(3)          ! q-vector
  REAL(DP) :: avg_iter
  INTEGER :: ig, i
  REAL(DP) :: xkold(3)
  REAL(DP), allocatable :: et_old(:,:)
  REAL(DP) :: rr, arg

  CALL start_clock( 'c_bands' )

  ! Initialize the diagonalization
  if (isolve == 1) then
    nbndx = nbnd ! CG
  elseif (isolve == 0) then
    nbndx = 4*nbnd ! Davidson
  else
    call errore('compute_u_kq', &
                'Don''t even try to use this isolve!', abs(isolve))
  endif
  cntrl_isolve = isolve
  max_cg_iter = 200
  iter = 1
  istep = 0
  ethr = conv_threshold
  lscf = .false.
  if (allocated(btype)) deallocate(btype)
  allocate(btype(nbndx,nkstot))
  btype(1:nbnd,:) = 1

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

  ! randomize a little bit
  do i = 1, nbnd
    do ig = 1, npw
      rr = 0.1d0*(2.d0*rndm() - 1.d0)
      arg = tpi * rndm()
      evc(ig,i) = evc(ig,i)*CMPLX(1.d0+rr*cos(arg),rr*sin(arg))
    enddo
  enddo

  ! Needed for LDA+U
  IF ( lda_plus_u ) CALL davcio( swfcatom, nwordatwfc, iunsat, ik, -1 )

  ! diagonalization of bands for k-point ik
  call diag_bands ( iter, ik, avg_iter )

!#ifdef __PARA
!  call mp_sum( avg_iter, inter_pool_comm )
!#endif
  avg_iter = avg_iter / nkstot

  !! debug
  !!write(stdout,'(5X,"ethr = ",1PE9.2,",  avg # of iterations =",0PF5.1)') &
  !!     ethr, avg_iter

  ! check if diagonalization was ok
  !!write(stdout,'(8F9.4)') et_old(1:nbnd,ik)*RytoeV
  !!write(stdout,'(8F9.4)') et(1:nbnd,ik)*RytoeV
  do i = 1, nbnd
    if (abs(et(i,ik) - et_old(i,ik))*RytoeV > 0.2d0) then
      write(stdout,'(5X,''ATTENTION: ik='',I4,''  ibnd='',I3,$)') ik, i
      write(stdout,'(2X,''eigenvalues differ too much!'')')
      write(stdout,'(5X,2(F10.4,2X))') et_old(i,ik)*RytoeV, et(i,ik)*RytoeV
    endif
  enddo

  ! restore the k-point and eigenvalues
  xk(:,ik) = xkold(:)
  et = et_old
  deallocate(et_old)

  ! restore wavefunctions
  evq = evc
  CALL get_buffer( evc, nwordwfc, iunwfc, ik)

  CALL stop_clock( 'c_bands' )  
  RETURN

END SUBROUTINE compute_u_kq
