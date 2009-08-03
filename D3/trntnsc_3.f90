!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine trntnsc_3 (phi, at, bg, iflg)
  !-----------------------------------------------------------------------
  !
  ! trasforms a COMPLEX third order tensor
  !(like the derivative of the dynamical matrix)
  ! from crystal to cartesian axis (iflg >=  1) or viceversa (iflg <= -1)
  !
  USE kinds, only : DP
  implicit none

  integer :: iflg
  ! input: gives the versus of the trans.

  complex (DP) :: phi (3, 3, 3)
  ! inp/out: the matrix to transform

  real (DP) :: at (3, 3), bg (3, 3)
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice

  integer :: i, j, k, l, m, n
  !
  !  counters on polarizations
  !
  complex (DP) :: wrk (3, 3, 3)
  ! a work array

  if (iflg.gt.0) then
     !
     ! forward transformation (crystal to cartesian axis)
     !

     call zcopy (27, phi, 1, wrk, 1)
     do m = 1, 3
        do i = 1, 3
           do j = 1, 3
              phi (m, i, j) = (0.d0, 0.d0)
              do n = 1, 3
                 do k = 1, 3
                    do l = 1, 3
                       phi (m, i, j) = phi (m, i, j) + wrk (n, k, l) * bg (i, k) &
                            * bg (j, l) * bg (m, n)
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
     do m = 1, 3
        do i = 1, 3
           do j = 1, 3
              wrk (m, i, j) = (0.d0, 0.d0)
              do n = 1, 3
                 do k = 1, 3
                    do l = 1, 3
                       wrk (m, i, j) = wrk (m, i, j) + phi (n, k, l) * at (k, i) &
                            * at (l, j) * at (n, m)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
     call zcopy (27, wrk, 1, phi, 1)
  endif
  return


end subroutine trntnsc_3
