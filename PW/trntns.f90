!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine trntns (phi, at, bg, iflg)
  !-----------------------------------------------------------------------
  !
  ! trasforms a tensor (like the stress)
  ! from crystal to cartesian axis (iflg >=  1) or viceversa (iflg <= -1)
  !
  USE kinds
  implicit none
  integer :: iflg, i, j, k, l
  real(DP) :: phi (3, 3), wrk (3, 3), at (3, 3), bg (3, 3)
  !
  if (iflg.gt.0) then
     !
     ! forward transformation (crystal to cartesian axis)
     !
     do i = 1, 3
        do j = 1, 3
           wrk (i, j) = phi (i, j)
        enddo
     enddo
     do i = 1, 3
        do j = 1, 3
           phi (i, j) = 0.0d0
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
           wrk (i, j) = 0.0d0
           do k = 1, 3
              do l = 1, 3
                 wrk (i, j) = wrk (i, j) + phi (k, l) * at (k, i) * at (l, j)
              enddo
           enddo
        enddo
     enddo
     do i = 1, 3
        do j = 1, 3
           phi (i, j) = wrk (i, j)
        enddo
     enddo
  endif
  !
  return
end subroutine trntns
