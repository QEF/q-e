!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine trntnsr_3 (phi, at, bg, iflg)
  !-----------------------------------------------------------------------
  !
  ! trasforms a REAL third order tensor 
  ! from crystal to cartesian axis (iflg >=  1) or viceversa (iflg <= -1)
  !
#include "f_defs.h"
  use kinds, only : DP
  implicit none

  integer :: iflg
  ! input: gives the versus of the trans.

  real(kind=DP) :: phi (3, 3, 3), at (3, 3), bg (3, 3)
  ! inp/out: the matrix to transform
  ! input: the direct lattice vectors 
  ! input: the reciprocal lattice

  integer :: i, j, k, l, m, n
  ! counters on polarizations
  real(kind=DP) :: wrk (3, 3, 3)
  ! a working array

  if(iflg.gt.0) then
     !
     ! forward transformation (crystal to cartesian axis)
     !
     call DCOPY(27, phi, 1, wrk, 1)

     do i = 1, 3
     do j = 1, 3
     do k = 1, 3
        phi (i, j, k) = 0.d0
        do l = 1, 3
        do m = 1, 3
        do n = 1, 3
           phi (i, j, k) = phi (i, j, k) +                     &
              wrk (l, m, n) * bg (i, l) * bg (j, m) * bg (k, n)
        enddo
        enddo
        enddo
     enddo
     enddo
     enddo
  else
     !
     ! backward transformation (cartesian to crystal axis)
     !
     do i = 1, 3
     do j = 1, 3
     do k = 1, 3
        wrk (i,j,k) = 0.d0
        do l = 1, 3
        do m = 1, 3
        do n = 1, 3
           wrk (i, j, k) = wrk (i, j, k) +                      &
              phi (l, m, n) * at (l, i) * at (m, j) * at (n, k)
        enddo
        enddo
        enddo
     enddo
     enddo
     enddo
     call DCOPY (27, wrk, 1, phi, 1)
  endif

  return
  end subroutine trntnsr_3
