!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine ccalbec (nkb, npwx, npw, nbnd, bec, vkb, psi)
  !-----------------------------------------------------------------------
  !
  !    This subroutine computes the dot product of the beta functions
  !    and the wavefunctions, and save them in the array bec.
  !
#include "machine.h"
  use parameters, only: DP
  implicit none
  !
  !   here the dummy variables
  !
  integer :: nkb, npwx, npw, nbnd
  ! input: the total number of beta functions
  ! input: the maximum number of plane waves
  ! input: the length of the vectors
  ! input: the number of bands
  complex(kind=DP) ::  vkb (npwx,nkb), psi (npwx,nbnd), bec (nkb,nbnd)
  ! input: the FT of the beta functions
  ! input: the wavefunctions
  ! output: dot product of the beta and the wavefunctions
  complex(kind=DP) :: alpha, beta
  !
  if (nkb.eq.0) return

  call start_clock ('ccalbec')
  alpha= (1.d0, 0.d0)
  beta = (0.d0, 0.d0)
  call ZGEMM ('C', 'N', nkb, nbnd, npw, alpha, vkb, npwx, psi, &
       npwx, beta, bec, nkb)
#ifdef PARA
  call reduce (2 * nkb * nbnd, bec)
#endif
  call stop_clock ('ccalbec')
  return
end subroutine ccalbec

