!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine trnvecc (u, at, bg, iflg)
  !-----------------------------------------------------------------------
  !
  !  transforms a COMPLEX vector in real space (like a displacement)
  !  from crystal to cartesian axis (iflag.gt.0) and viceversa (iflag.le.0
  !
  USE kinds, only : DP
  implicit none

  integer :: iflg
  ! input: gives the versus of the transformatio

  real(DP) :: at (3, 3), bg (3, 3)
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice vectors

  complex(DP) :: u (3)
  ! inp/out: the vector to transform

  integer :: i, k
  !   !  counter on polarizations

  !/

  complex(DP) :: wrk (3)
  ! auxiliary variable
  if (iflg.gt.0) then
     !
     ! forward transformation :
     !
     do i = 1, 3
        wrk (i) = u (i)

     enddo
     do i = 1, 3
        u (i) = 0.d0
        do k = 1, 3
           u (i) = u (i) + wrk (k) * at (i, k)
        enddo
     enddo
  else
     !
     ! backward transformation :
     !
     do i = 1, 3
        wrk (i) = 0.d0
        do k = 1, 3
           wrk (i) = wrk (i) + u (k) * bg (k, i)
        enddo

     enddo
     do i = 1, 3
        u (i) = wrk (i)
     enddo
  endif
  return
end subroutine trnvecc
