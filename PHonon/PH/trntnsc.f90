!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine trntnsc (phi, at, bg, iflg)
  !-----------------------------------------------------------------------
  !
  ! trasforms a COMPLEX tensor (like the dynamical matrix)
  ! from crystal to cartesian axis (iflg >=  1) or viceversa (iflg <= -1)
  !
  USE kinds, only : DP
  implicit none

  integer :: iflg
  ! input: gives the versus of the trans.

  complex(DP) :: phi (3, 3)
  ! inp/out: the matrix to transform

  real(DP) :: at (3, 3), bg (3, 3)
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice

  integer :: i, j, k, l
  !
  !  counters on polarizations
  ! /
  !/


  complex(DP) :: wrk (3, 3)
  ! a working array
  if (iflg.gt.0) then
     !
     ! forward transformation (crystal to cartesian axis)
     !

     call zcopy (9, phi, 1, wrk, 1)
     do i = 1, 3
        do j = 1, 3
           phi (i, j) = (0.d0, 0.d0)
           do k = 1, 3
              do l = 1, 3
                 phi (i, j) = phi (i, j) + wrk (k, l) * bg (i, k) * bg (j, l)
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
           wrk (i, j) = (0.d0, 0.d0)
           do k = 1, 3
              do l = 1, 3
                 wrk (i, j) = wrk (i, j) + phi (k, l) * at (k, i) * at (l, j)
              enddo
           enddo
        enddo
     enddo
     call zcopy (9, wrk, 1, phi, 1)
  endif
  return
end subroutine trntnsc
