!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
function dsum (n, vect, inc)
  !-----------------------------------------------------------------------
  !
  !     This function compute the sum of all the elements of a vector vect
  !
  USE kinds
  implicit none
  !
  !    first the dummy variables
  !
  integer :: n, inc
  ! input: dimension of the vector
  ! input: distance between the elements
  real(kind=DP) :: vect (n), dsum
  ! input: the vector
  ! output: the sum of the elements
  !
  !   local variables
  !


  integer :: i
  ! counter on the elements
  dsum = 0.d0
  if (n.lt.0.or.inc.le.0) return
  do i = 1, n, inc
     dsum = dsum + vect (i)

  enddo
  return
end function dsum
