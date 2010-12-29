!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine bachel (alps, aps, npseu, lmax)
  !----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants , ONLY : pi
  implicit none
  !
  !   First I/O variables
  !
  integer :: npseu, lmax (npseu)
  ! input:  number of pseudopotential
  ! input:  max. angul. momentum of the ps
  real(DP) :: alps (3, 0:3, npseu), aps (6, 0:3, npseu)
  ! input:  the b_l coefficient
  ! in/out: the a_l coefficient
  !
  !   Here local variables
  !
  integer :: np, lmx, l, i, j, k, ia, ka, nik
  ! counter on number of pseudopot.
  ! aux. var. (max. ang. mom. of a fix. ps
  ! counter on angular momentum

  real(DP) :: s (6, 6), alpl, alpi, ail
  ! auxiliary array
  ! first real aux. var. (fix. value of al
  ! second real aux. var. (fix. value of a
  ! third real aux. var.
  !
  do np = 1, npseu
     lmx = lmax (np)
     do l = 0, lmx
        do k = 1, 6
           ka = mod (k - 1, 3) + 1
           alpl = alps (ka, l, np)
           do i = 1, k
              ia = mod (i - 1, 3) + 1
              alpi = alps (ia, l, np)
              ail = alpi + alpl
              s (i, k) = sqrt (pi / ail) / 4.d0 / ail
              nik = int ( (k - 1) / 3) + int ( (i - 1) / 3) + 1
              do j = 2, nik
                 s (i, k) = s (i, k) / 2.d0 / ail * (2 * j - 1)
              enddo
           enddo
        enddo
        !
        do i = 1, 6
           do j = i, 6
              do k = 1, i - 1
                 s (i, j) = s (i, j) - s (k, i) * s (k, j)
              enddo
              if (i.eq.j) then
                 s (i, i) = sqrt (s (i, i) )
              else
                 s (i, j) = s (i, j) / s (i, i)
              endif
           enddo
        enddo
        !
        aps (6, l, np) = - aps (6, l, np) / s (6, 6)
        do i = 5, 1, - 1
           aps (i, l, np) = - aps (i, l, np)
           do k = i + 1, 6
              aps (i, l, np) = aps (i, l, np) - aps (k, l, np) * s (i, k)
           enddo
           aps (i, l, np) = aps (i, l, np) / s (i, i)
        enddo
     enddo

  enddo
  return
end subroutine bachel
