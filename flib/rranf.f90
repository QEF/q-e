!
! Copyright (C) 2002-2003 CP group
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

