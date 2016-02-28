!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine setqmod (ngm, xq, g, qmod, qpg)
  !-----------------------------------------------------------------------
  !
  ! This subroutine puts in qmod the modulus of q+G for the interpolation
  ! table used to compute qgm
  !
  USE kinds, only : DP
  
  implicit none

  integer :: ngm
  ! input: the number of G vectors

  real(DP) :: xq (3), g (3, ngm), qmod (ngm), qpg (3, ngm)
  ! input: the q vector
  ! input: the G vectors
  ! output: the modulus of the G vectors
  ! output: the q+G vectors

  integer :: ig
  ! counter on G vectors
  
  do ig = 1, ngm
     qpg (1, ig) = xq (1) + g (1, ig)
     qpg (2, ig) = xq (2) + g (2, ig)
     qpg (3, ig) = xq (3) + g (3, ig)
     qmod (ig) = qpg (1, ig) **2 + qpg (2, ig) **2 + qpg (3, ig) **2
  enddo

  return

end subroutine setqmod
