!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------

subroutine recips (a1, a2, a3, b1, b2, b3)
  !---------------------------------------------------------------------
  !! This routine generates the reciprocal lattice vectors \(b_1,b_2,b_3\)
  !! given the real space vectors \(a_1,a_2,a_3\). The \(b\)'s are units of
  !! \(2 \pi/a\).
  !
  use kinds, ONLY: DP
  implicit none
  real(DP) :: a1 (3)
  !! input: first direct lattice vector
  real(DP) :: a2 (3)
  !! input: second direct lattice vector
  real(DP) :: a3 (3)
  !! input: third direct lattice vector
  real(DP) :: b1 (3)
  !! output: first reciprocal lattice vector
  real(DP) :: b2 (3)
  !! output: second reciprocal lattice vector
  real(DP) :: b3 (3)
  !! output: third reciprocal lattice vector
  !
  ! ... local variables
  !
  real(DP) :: den, s
  ! the denominator
  ! the sign of the permutations
  integer :: iperm, i, j, k, l, ipol
  ! counter on the permutations
  !\
  !  Auxiliary variables
  !/
  !
  ! Counter on the polarizations
  !
  !    first we compute the denominator
  !
  den = 0
  i = 1
  j = 2
  k = 3
  s = 1.d0
100 do iperm = 1, 3
     den = den + s * a1 (i) * a2 (j) * a3 (k)
     l = i
     i = j
     j = k
     k = l
  enddo
  i = 2
  j = 1
  k = 3
  s = - s
  if (s.lt.0.d0) goto 100
  !
  !    here we compute the reciprocal vectors
  !
  i = 1
  j = 2
  k = 3
  do ipol = 1, 3
     b1 (ipol) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) ) / den
     b2 (ipol) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) ) / den
     b3 (ipol) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) ) / den
     l = i
     i = j
     j = k
     k = l
  enddo
  return
end subroutine recips
