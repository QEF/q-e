! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE greenfunction(ik, psi, g_psi, q)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the Green function operator
  ! ... (formula here)
  ! ... We use Hartree atomic units; since G is the inverse of an
  ! ... energy: G => G / ryd_to_hartree
  !
  USE kinds,                       ONLY : DP
  USE io_global,                   ONLY : stdout  
  USE becmod,                      ONLY : becp
  USE wavefunctions_module,        ONLY : evc
  USE noncollin_module,            ONLY : npol
  USE pwcom
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE gipaw_module

  !-- parameters ---------------------------------------------------------
  IMPLICIT none
  INTEGER, INTENT(IN) :: ik
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx,nbnd)  ! psi is changed on output!!!
  COMPLEX(DP), INTENT(OUT) :: g_psi(npwx,nbnd)
  REAL(DP) :: q(3)

  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  complex(dp), allocatable :: ps(:,:), work (:)
  real(dp), allocatable :: h_diag (:,:), eprec (:)
  real(dp) :: anorm, thresh, gk(3), dxk(3)
  integer :: ibnd, ig, lter
  logical :: conv_root, q_is_zero
  complex(dp), external :: ZDOTC
  external ch_psi_all, cg_psi
 
  ! start clock
  call start_clock ('greenf')

  ! allocate memory
  allocate (work(npwx), ps(nbnd,nbnd), h_diag(npwx,nbnd), &
            eprec(nbnd), becp(nkb,nbnd))

  ! check if |q| is zero
  q_is_zero = .false.
  if (sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3)) < 1d-8) q_is_zero = .true.  
  if (q_is_zero) evq(:,:) = evc(:,:)

  !====================================================================
  ! apply -Q_{k+q} to the r.h.s.
  !====================================================================
  ! project on <evq|: ps(i,j) = <evq(i)|psi(j)>
  ps = (0.d0,0.d0)
  CALL ZGEMM('C', 'N', nbnd_occ (ik), nbnd_occ (ik), npw, &
             (1.d0,0.d0), evq(1,1), npwx, psi(1,1), npwx, (0.d0,0.d0), &
             ps(1,1), nbnd)
#ifdef __PARA
  call reduce (2 * nbnd * nbnd_occ (ik), ps(1,1))
#endif

  !! this is the case with overlap (ultrasoft)
  !! g_psi is used as work space to store S|evc>
  !!
  !!CALL ccalbec (nkb, npwx, npw, nbnd_occ(ik), becp, vkb, evc)
  !!CALL s_psi (npwx, npw, nbnd_occ(ik), evc, g_psi)
  !! |psi> = -(|psi> - S|evc><evc|psi>)
  !!
  !!CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
  !!     (1.d0,0.d0), g_psi(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
  !!     psi(1,1), npwx )

  ! |psi> = -(1 - |evq><evq|) |psi>
  CALL ZGEMM('N', 'N', npw, nbnd_occ(ik), nbnd_occ(ik), &
             (1.d0,0.d0), evq(1,1), npwx, ps(1,1), nbnd, (-1.d0,0.d0), &
             psi(1,1), npwx)


  !====================================================================
  ! solve the linear system (apply G_{k+q})
  !====================================================================
  ! convergence treshold
  thresh = sqrt(conv_threshold)   ! sqrt(of that of PARATEC)

  ! use the hamiltonian at k+q
  do ig = 1, npw
    gk(1) = (xk(1,ik) + g(1,igk(ig)) + q(1)) * tpiba
    gk(2) = (xk(2,ik) + g(2,igk(ig)) + q(2)) * tpiba
    gk(3) = (xk(3,ik) + g(3,igk(ig)) + q(3)) * tpiba
    g2kin (ig) = gk(1)**2 + gk(2)**2 + gk(3)**2
  enddo

  ! preconditioning of the linear system
  do ibnd = 1, nbnd_occ (ik)
     do ig = 1, npw
        work (ig) = g2kin (ig) * evq (ig, ibnd)
     enddo
     eprec (ibnd) = 1.35d0 * ZDOTC (npw, evq (1, ibnd), 1, work, 1)
  enddo
#ifdef __PARA
  call reduce (nbnd_occ (ik), eprec)
#endif
  do ibnd = 1, nbnd_occ (ik)
     do ig = 1, npw
        h_diag (ig, ibnd) = 1.d0 / max (1.0d0, g2kin (ig) / eprec (ibnd) )
     enddo
  enddo

  if (.not. q_is_zero) then
    dxk = xk(:,ik) + q
    call init_us_2(npw, igk, dxk, vkb)
  else
    call init_us_2(npw, igk, xk(1,ik), vkb)
  endif
  call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, psi)
    
  ! initial guess
  g_psi(:,:) = (0.d0, 0.d0)

  ! solve linear system  
  conv_root = .true.
  call cgsolve_all (ch_psi_all, cg_psi, et(1,ik), psi, g_psi, &
       h_diag, npwx, npw, thresh, ik, lter, conv_root, anorm, &
       nbnd_occ(ik) )

  !! debug  
  !!write(stdout, '(5X,''cgsolve_all converged in '',I3,'' iterations'')') &
  !!      lter

  if (.not.conv_root) WRITE( stdout, '(5x,"ik",i4," ibnd",i4, &
       & " linter: root not converged ",e10.3)') &
       ik, ibnd, anorm

  ! convert to Hartree
  g_psi(:,:) = g_psi(:,:) / ryd_to_hartree

  call flush_unit( stdout )
  call stop_clock('greenf')
 
  ! free memory
  deallocate (work, h_diag, eprec, ps, becp)

END SUBROUTINE greenfunction
