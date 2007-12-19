!
! Copyright (C) 2003-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine H_h(e,h,Ah)
  !-----------------------------------------------------------------------
#include "f_defs.h"
  !
  USE kinds, only: DP
  USE wvfct, ONLY: nbnd, npwx, npw, g2kin, igk
  USE gvect, ONLY : gstart
  USE uspp,     ONLY : vkb, nkb
  USE lsda_mod, ONLY : current_spin
  USE scf,      ONLY : vrs
  use becmod, only: rbecp, calbec
  use cgcom
  !
  implicit none
  !
  real(DP):: e(nbnd)
  complex(DP):: h(npwx,nbnd), Ah(npwx,nbnd)
  !
  integer:: j,ibnd
  !
  call start_clock('h_h')
  !
  ! [(k+G)^2 - e ]psi
  do ibnd = 1,nbnd
     ! set to zero the imaginary part of h at G=0
     !  needed for numerical stability
     if (gstart==2) h(1,ibnd) = CMPLX( DBLE(h(1,ibnd)),0.d0)
     do j = 1,npw
        ah(j,ibnd) = (g2kin(j)-e(ibnd)) * h(j,ibnd)
     end do
  end do
  ! V_Loc psi
  call vloc_psi(npwx, npw, nbnd, h, vrs(1,current_spin), ah)
  ! V_NL psi
   call calbec  ( npw, vkb, h, rbecp )
  if (nkb > 0) call add_vuspsi (npwx, npw, nbnd, h, ah)
  ! set to zero the imaginary part of ah at G=0
  !  needed for numerical stability
  if (gstart==2) then
     do ibnd = 1, nbnd
        ah(1,ibnd) = CMPLX( DBLE(ah(1,ibnd)),0.d0)
     end do
  end if
  !
  call stop_clock('h_h')
  !
  return
end subroutine H_h
