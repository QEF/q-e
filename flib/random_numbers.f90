!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE random_numbers
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, SAVE :: irand
  !
  INTERFACE gauss_dist
     !
     MODULE PROCEDURE gauss_dist_scal, gauss_dist_vect
     !
  END INTERFACE gauss_dist
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    FUNCTION rranf()
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP) :: rranf
      !
      INTEGER :: kk, m, konst
      DATA m/100001/, konst/125/
      SAVE m, konst
      !
      m = m * konst
      m = m - 2796203 * ( m / 2796203 )
      !
      rranf = DBLE( m ) / 2796203.D0
      !
      RETURN
      !
    END FUNCTION rranf
    !
    !------------------------------------------------------------------------
    FUNCTION randy()
      !------------------------------------------------------------------------
      !
      ! ... use machine-specific random-number generator when available
      !
      REAL(DP) :: randy
#ifdef __AIX
#define __USE_SYSTEM_RAND
      randy = rand()
#endif
      !
      ! ... use fortran random-number generator in all other cases
      !
#ifndef __USE_SYSTEM_RAND
      INTEGER , PARAMETER  :: m    = 714025, &
                              ia   = 1366, &
                              ic   = 150889, &
                              ntab = 97
      REAL(DP), PARAMETER  :: rm = 1.D0 / m
      INTEGER              :: ir(ntab), iff, idum, j, iy
      DATA iff /0/, idum/0/
      SAVE iff, idum, iy, ir
      !
      !
      IF ( idum < 0 .OR. iff == 0 ) THEN
         !
         iff = 1
         idum = MOD( ic - idum, m )
         !
        DO j=1,ntab
           idum=mod(ia*idum+ic,m)
           ir(j)=idum
        END DO
        idum=mod(ia*idum+ic,m)
        iy=idum
      END IF
      j=1+(ntab*iy)/m
      IF(j.gt.ntab.or.j.lt.1) call errore('randy','j out of range',j)
      iy=ir(j)
      randy=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
#endif
      RETURN
      !
    END FUNCTION randy
    !
    !------------------------------------------------------------------------
    function rndm()
      !------------------------------------------------------------------------
      !
      ! ... random number generator equivalent to ran1 of Num.Rec.
      !
      IMPLICIT NONE
      !
      REAL(DP) :: rndm
      REAL(DP) :: shuffle(32)
      INTEGER  :: i
      LOGICAL  :: first
      DATA first / .TRUE. /
      SAVE first, shuffle, i
      !
      ! ... starting seed, must be not be 0
      !
      IF ( first ) irand = -1
      !
      IF ( first .OR. irand < 0 ) THEN
         !
         irand = - irand
         !
         DO i = 32 + 8, 1, - 1
            !
            shuffle( MIN( i, 32 ) ) = rndx( irand )
            !
         END DO
         !
         i = 32 * shuffle(1) + 1
         !
         first = .FALSE.
         !
      end if
      !
      rndm = shuffle(i)
      !
      shuffle(i) = rndx( irand )
      !
      i = 32 * rndm + 1
      !
      RETURN
      !
    END FUNCTION rndm
    !
    !------------------------------------------------------------------------
    FUNCTION rndx( irand )
      !------------------------------------------------------------------------
      !
      ! ... random number generator equivalent to ran0 of Num.Rec.
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(INOUT) :: irand
      REAL(DP)               :: rndx
      !
      INTEGER  :: im, ia, iq, ir, is, it
      REAL(DP) :: obm
      LOGICAL  :: first
      DATA first / .TRUE. /
      SAVE im, ia, iq, ir, obm, first
      !
      IF ( first ) THEN
         !
         ! ... this is 2**31-1 avoiding overflow
         !
         im  = 2 * ( 2**30 - 1 ) + 1
         obm = 1.0 / im
         ia  = 7*7*7*7*7
         iq  = im / ia
         ir  = im - ia * iq
         !
         first = .FALSE.
         !
      END IF
      !
      is = irand / iq
      it = irand - is * iq
      !
      irand = ia * it - is * ir
      !
      IF ( irand < 0 ) irand = irand + im
      !
      rndx = irand * obm
      !
      RETURN
      !
    END FUNCTION rndx
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_rndm_seed( iseed )
      !------------------------------------------------------------------------
      !
      ! ... this subroutine initialize the random number with the given seed
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: iseed
      INTEGER             :: irand
      !
      REAL(DP) :: dummy
      !
      IF ( iseed < 0 ) &
         CALL errore( 'set_rndm_seed', 'seed should be a positive integer', 1 )
      !
      ! ... make sure rndm() has been called already once !
      !
      dummy = rndm()
      irand = - iseed
      !
      RETURN
      !
    END SUBROUTINE set_rndm_seed
    !
    !-----------------------------------------------------------------------
    FUNCTION gauss_dist_scal( mu, sigma )
      !-----------------------------------------------------------------------
      !
      ! ... this function generates a number taken from a normal
      ! ... distribution of mean value \mu and variance \sigma
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: mu
      REAL(DP), INTENT(IN) :: sigma
      REAL(DP)             :: gauss_dist_scal
      !
      REAL(DP) :: x1, x2, w, coeff
      !
      !
      gaussian_loop: DO
         !
         x1 = 2.D0 * rndm() - 1.D0
         x2 = 2.D0 * rndm() - 1.D0
         !
         w = x1 * x1 + x2 * x2
         !
         IF ( w < 1.D0 ) EXIT gaussian_loop
         !
      END DO gaussian_loop
      !
      w = SQRT( ( - 2.D0 * LOG( w ) ) / w )
      !
      gauss_dist_scal = x1 * w * sigma + mu
      !
      RETURN
      !
    END FUNCTION gauss_dist_scal
    !    
    !-----------------------------------------------------------------------
    FUNCTION gauss_dist_vect( mu, sigma, dim )
      !-----------------------------------------------------------------------
      !
      ! ... this function generates an array of numbers taken from a normal
      ! ... distribution of mean value \mu and variance \sigma
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: mu
      REAL(DP), INTENT(IN) :: sigma
      INTEGER,  INTENT(IN) :: dim
      REAL(DP)             :: gauss_dist_vect( dim )
      !
      REAL(DP) :: x1, x2, w, coeff
      INTEGER  :: i
      !
      !
      DO i = 1, dim, 2
         !
         gaussian_loop: DO
            !
            x1 = 2.D0 * rndm() - 1.D0
            x2 = 2.D0 * rndm() - 1.D0
            !
            w = x1 * x1 + x2 * x2
            !
            IF ( w < 1.D0 ) EXIT gaussian_loop
            !
         END DO gaussian_loop
         !
         w = SQRT( ( - 2.D0 * LOG( w ) ) / w )
         !
         gauss_dist_vect(i) = x1 * w * sigma
         !
         IF ( i >= dim ) EXIT
         !
         gauss_dist_vect(i+1) = x2 * w * sigma
         !
      END DO
      !
      gauss_dist_vect(:) = gauss_dist_vect(:) + mu
      !
      RETURN
      !
    END FUNCTION gauss_dist_vect
    !
END MODULE random_numbers
