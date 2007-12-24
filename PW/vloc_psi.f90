!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine vloc_psi(lda, n, m, psi, v, hpsi)
  !-----------------------------------------------------------------------
  !
  USE kinds, only : DP
  USE gsmooth, ONLY : nls, nlsm, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE wvfct,   ONLY : igk
  USE wavefunctions_module,  ONLY: psic
  USE control_flags, ONLY : use_task_groups
  USE mp_global,     ONLY : nogrp, ogrp_comm, me_pool
  USE task_groups,   ONLY : nolist
  !
  implicit none
  !
  integer :: lda, n, m
  complex(DP) :: psi (lda, m), hpsi (lda, m)
  real(DP) :: v(nrxxs)
  !
  complex(DP) :: fp, fm
  integer :: ibnd, j, incr
  ! counters

  call start_clock ('vloc_psi')
  !
  IF( use_task_groups ) THEN
     incr = 2 * nogrp
  ELSE
     incr = 2
  END IF
  !
  ! the local potential V_Loc psi. First bring psi to real space
  !
  do ibnd = 1, m, incr
     psic(:) = (0.d0, 0.d0)
     if (ibnd < m) then
        ! two ffts at the same time
        do j = 1, n
           psic (nls (igk(j))) =       psi(j, ibnd) + (0.0d0,1.d0)*psi(j, ibnd+1)
           psic (nlsm(igk(j))) = CONJG(psi(j, ibnd) - (0.0d0,1.d0)*psi(j, ibnd+1))
        enddo
     else
        do j = 1, n
           psic (nls (igk(j))) =       psi(j, ibnd)
           psic (nlsm(igk(j))) = CONJG(psi(j, ibnd))
        enddo
     end if
     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
     !
     !   product with the potential v on the smooth grid
     !
     do j = 1, nrxxs
        psic (j) = psic (j) * v(j)
     enddo
     !
     !   back to reciprocal space
     !
     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
     !
     !   addition to the total product
     !
     if (ibnd < m) then
        ! two ffts at the same time
        do j = 1, n
           fp = (psic (nls(igk(j))) + psic (nlsm(igk(j))))*0.5d0
           fm = (psic (nls(igk(j))) - psic (nlsm(igk(j))))*0.5d0
           hpsi (j, ibnd)   = hpsi (j, ibnd)   + CMPLX( DBLE(fp), AIMAG(fm))
           hpsi (j, ibnd+1) = hpsi (j, ibnd+1) + CMPLX(AIMAG(fp),- DBLE(fm))
        enddo
     else
        do j = 1, n
           hpsi (j, ibnd)   = hpsi (j, ibnd)   + psic (nls(igk(j)))
        enddo
     end if
  enddo
  call stop_clock ('vloc_psi')
  return
end subroutine vloc_psi
