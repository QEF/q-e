!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
function rndm ()
  !
  !   RANDOM NUMBER GENERATOR equivalent to ran1 of Num.Rec.
  !
  use parameters
  implicit none
  integer :: irand
  common/random_number/ irand
  real(kind=DP) :: rndm, rnd, shuffle (32)
  integer :: i
  logical :: first
  data first / .true. /
  save first, shuffle, i
  if (first.or.irand.lt.0) then
     irand = - irand
     if (first) irand=1 ! starting seed, must be not be 0
     do i = 32 + 8, 1, - 1
        shuffle (min (i, 32) ) = rnd(irand)
     enddo
     i = 32 * shuffle (1) + 1
     first = .false.
  endif
  rndm = shuffle (i)
  shuffle (i) = rnd(irand)

  i = 32 * rndm + 1
  return

end function rndm

function rnd (irand)
  !
  !   RANDOM NUMBER GENERATOR equivalent to ran0 of Num.Rec.
  !
  use parameters
  implicit none
  integer :: im, ia, iq, ir, irand, is, it
  real(kind=DP) :: rnd, obm
  logical :: first
  data first / .true. /
  save im, ia, iq, ir, obm, first
  if (first) then
     ! this is 2**31-1 avoiding overflow
     im = 2 * (2**30 - 1) + 1
     obm = 1.0 / im
     ia = 7**5
     iq = im / ia
     ir = im - ia * iq
     ! starting seed, must be not be 0
!     irand = 1
     first = .false.
  endif
  is = irand / iq
  it = irand-is * iq
  irand = ia * it - is * ir
  if (irand.lt.0) irand = irand+im
  rnd = irand * obm
  return
end function rnd


subroutine set_rndm_seed(iseed)
!
! this subroutine initialize the random number with the given seed
!
  use parameters
  implicit none
  integer :: irand,iseed
  common/random_number/ irand
  real(kind=DP) :: dummy, rndm
  if (iseed.le.0) call error('set_rndm_seed', &
                  'seed should be a positive integer',1)
! make sure rndm() has been called already once !
  dummy = rndm()
  irand = - iseed

  return
end subroutine set_rndm_seed

