!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine swap (n, x, x1)
  !-----------------------------------------------------------------------
  ! swap array x with array x1
  !
  USE kinds
  implicit none
  !
  !   I/O variables
  !
  integer :: n
  ! input: dimension of the vector
  real(DP) :: x (n), x1 (n)
  ! I/O: the vectors
  !
  !   local variables
  !
  integer :: i
  ! counter on the elements

  real(DP) :: xswap
  ! work
  do i = 1, n
     xswap = x (i)
     x (i) = x1 (i)
     x1 (i) = xswap

  enddo
  return
end subroutine swap
