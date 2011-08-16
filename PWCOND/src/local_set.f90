!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine local_set(n1, n2, n3, n4, n5, n6, n7, n8)
!
! To set up the number of all the orbitals to be 0 for
! local potential calculations
!
  implicit none
  integer :: n1, n2, n3, n4, n5, n6, n7, n8

  n1 = 0
  n2 = 0
  n3 = 0
  n4 = 0
  n5 = 0
  n6 = 0
  n7 = 0
  n8 = 0

  return
end subroutine local_set

