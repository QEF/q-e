!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine drhodv(nu_i)
  !-----------------------------------------------------------------------
  !
  !  calculate the electronic term <psi|dv|dpsi> of the dynamical matrix
  !
#include "f_defs.h"
  use pwcom
  use cgcom
  implicit none
  integer :: nu_i
  !
  integer :: nu_j, ibnd, kpoint
  real(DP) :: dynel(nmodes), work(nbnd)
  !
  call start_clock('drhodv')
  !
  dynel(:) = 0.d0
  kpoint = 1
  ! do kpoint=1,nks
  !
  !** calculate the dynamical matrix (<DeltaV*psi(ion)|\DeltaPsi(ion)>)
  !
  do nu_j = 1,nmodes
     !
     ! DeltaV*psi(ion) for mode nu_j is recalculated
     !
     call dvpsi_kb(kpoint,nu_j)
     !
     !     this is the real part of <DeltaV*Psi(ion)|DeltaPsi(ion)>
     !
     call pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi ,npwx,work)
     do ibnd = 1,nbnd
        dynel(nu_j) = dynel(nu_j) + 2.0d0*wk(kpoint)*work(ibnd)
     end do
  end do
#ifdef __PARA
  call reduce(nmodes,dynel)
#endif
  !
  ! NB this must be done only at the end of the calculation!
  !
  do nu_j = 1,nmodes
     dyn(nu_i,nu_j) = - (dyn(nu_i,nu_j)+dynel(nu_j))
  end do
  !
  call stop_clock('drhodv')
  !
  return
end subroutine drhodv

