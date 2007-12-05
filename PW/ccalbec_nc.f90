!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine ccalbec_nc (nkb, npwx, npw, npol, nbnd, bec, vkb, psi)
  !-----------------------------------------------------------------------
  !
  !    This subroutine computes the dot product of the beta functions
  !    and the wavefunctions, and save them in the array bec.
  !
#include "f_defs.h"
  USE kinds, only: DP
  implicit none
  !
  !   here the dummy variables
  !
  integer :: nkb, npwx, npw, nbnd, npol
  ! input: the total number of beta functions
  ! input: the maximum number of plane waves
  ! input: the length of the vectors
  ! input: the number of bands
  ! input: number of coorrdinates of wfc
  complex(DP) :: vkb (npwx,nkb),psi(npwx,npol,nbnd),bec(nkb,npol,nbnd)
  ! input: the FT of the beta functions
  ! input: the wavefunctions
  ! output: dot product of the beta and the wavefunctions
  complex(DP) :: alpha, beta
  complex(DP), external :: ZDOTC
  !
  if (nkb.eq.0) return

  call start_clock ('calbec')
  alpha= (1.d0, 0.d0)
  beta = (0.d0, 0.d0)
  call ZGEMM ('C', 'N', nkb, nbnd*npol, npw, alpha, vkb, npwx, psi, &
       npwx, beta, bec, nkb)

#ifdef __PARA
  call reduce (2 * nkb * nbnd*npol, bec)
#endif
  call stop_clock ('calbec')
  return
end subroutine ccalbec_nc
