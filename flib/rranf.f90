!
! Copyright (C) 2002-2004 PWSCF-CP-FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
!
! random number generation
! ========================
!
!-------------------------------------------------------------------------

      FUNCTION RRANF()
        USE kinds
        IMPLICIT NONE
        REAL(dbl) :: RRANF
        INTEGER :: KK, M, KONST
        DATA M/100001/, KONST/125/
        SAVE M, KONST
          M=M*KONST
          M=M-2796203 * (M/2796203)
          RRANF = REAL(M)/2796203.D0
        RETURN
      END FUNCTION

!-------------------------------------------------------------------------
      real(kind=8) function randy()
!-------------------------------------------------------------------------
!
! Use machine-specific random-number generator when available
!
#ifdef __CRAYY
#define __USE_SYSTEM_RAND
      randy = ranf()
#endif
#ifdef __SX4
#define __USE_SYSTEM_RAND
      randy=random(0)
#endif
#ifdef __AIX
#define __USE_SYSTEM_RAND
      randy=rand()
#endif
!
! Use fortran random-number generator in all other cases
!
#ifndef __USE_SYSTEM_RAND
      integer m, ia, ic, ntab
      real(kind=8) rm
      parameter (ntab=97,m=714025,ia=1366,ic=150889,rm=1.0/m)
      integer ir(ntab), iff, idum, j, iy
      data iff /0/, idum/0/
      save iff, idum, iy, ir
!
!
      if(idum.lt.0.or.iff.eq.0) then
        iff=1
        idum=mod(ic-idum,m)
        do j=1,ntab
           idum=mod(ia*idum+ic,m)
           ir(j)=idum
        end do
        idum=mod(ia*idum+ic,m)
        iy=idum
      endif
      j=1+(ntab*iy)/m
      if(j.gt.ntab.or.j.lt.1) call errore('randy','j out of range',j)
      iy=ir(j)
      randy=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
#endif
      return
      end
!
function rndm ()
  !
  !   RANDOM NUMBER GENERATOR equivalent to ran1 of Num.Rec.
  !
  USE kinds
  implicit none
  integer :: irand
  common/random_number/ irand
  real(kind=DP) :: rndm, shuffle (32)
  real(kind=DP), external :: rndx
  integer :: i
  logical :: first
  data first / .true. /
  save first, shuffle, i
  !
  if (first) irand = -1 ! starting seed, must be not be 0
  !
  if (first.or.irand.lt.0) then
     irand = - irand
     do i = 32 + 8, 1, - 1
        shuffle (min (i, 32) ) = rndx (irand)
     enddo
     i = 32 * shuffle (1) + 1
     first = .false.
  endif
  rndm = shuffle (i)
  shuffle (i) = rndx (irand)

  i = 32 * rndm + 1
  return
 
end function rndm

function rndx (irand)
  !
  !   RANDOM NUMBER GENERATOR equivalent to ran0 of Num.Rec.
  !
  USE kinds
  implicit none
  integer :: im, ia, iq, ir, irand, is, it
  real(kind=DP) :: rndx, obm
  logical :: first
  data first / .true. /
  save im, ia, iq, ir, obm, first
  !
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
  rndx = irand * obm
  return
end function rndx


subroutine set_rndm_seed(iseed)
!
! this subroutine initialize the random number with the given seed
!
  USE kinds
  implicit none
  integer :: irand,iseed
  common/random_number/ irand
  real(kind=DP) :: dummy, rndm
  if (iseed.le.0) call errore('set_rndm_seed', &
                  'seed should be a positive integer',1)
! make sure rndm() has been called already once !
  dummy = rndm()
  irand = - iseed

  return
end subroutine set_rndm_seed

