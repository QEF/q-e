!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine H_h(e,h,Ah)
  !-----------------------------------------------------------------------
#include "machine.h"
  !
  use parameters, only: DP
  use pwcom
  use rbecmod
  use cgcom
  !
  implicit none
  !
  real(kind=DP):: e(nbnd)
  complex(kind=DP):: h(npwx,nbnd), Ah(npwx,nbnd)
  !
  integer:: j,ibnd
  !
  call start_clock('h_h')
  !
  ! [(k+G)^2 - e ]psi
  do ibnd = 1,nbnd
     ! set to zero the imaginary part of h at G=0
     !  needed for numerical stability
     if (gstart==2) h(1,ibnd) = cmplx(DREAl(h(1,ibnd)),0.d0)
     do j = 1,npw
        ah(j,ibnd) = (g2kin(j)-e(ibnd)) * h(j,ibnd)
     end do
  end do
  ! V_Loc psi
  call vloc_psi(npwx, npw, nbnd, h, vrs(1,current_spin), ah)
  ! V_NL psi
   call pw_gemm ('Y', nkb, nbnd, npw, vkb, npwx, h, npwx, becp, nkb)
  if (nkb.gt.0) call add_vuspsi (npwx, npw, nbnd, h, ah)
  ! set to zero the imaginary part of ah at G=0
  !  needed for numerical stability
  if (gstart==2) then
     do ibnd = 1, nbnd
        ah(1,ibnd) = cmplx(DREAL(ah(1,ibnd)),0.d0)
     end do
  end if
  !
  call stop_clock('h_h')
  !
  return
end subroutine H_h
