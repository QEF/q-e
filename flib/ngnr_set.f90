!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

   SUBROUTINE ngnr_set( alat, a1, a2, a3, gcut, qk, ng, nr1, nr2, nr3 )

!  this routine calculates the storage required for G vectors arrays
!  ----------------------------------------------
!  END manual

! ... declare modules
      USE kinds, ONLY: dbl
      USE mp, ONLY: mp_max, mp_min, mp_sum
      USE mp_global, ONLY: mpime, nproc, group

      IMPLICIT NONE

      INTEGER, INTENT(OUT) :: nr1, nr2, nr3, ng
      REAL(dbl), INTENT(IN) :: alat, a1(3), a2(3), a3(3), gcut, qk(3)

! ... declare other variables
      INTEGER :: i, j, k
      INTEGER :: nr1x, nr2x, nr3x
      INTEGER :: nr1tab, nr2tab, nr3tab, nr
      INTEGER :: nb(3)
      REAL(dbl) :: gsq, sqgc
      REAL(dbl) :: c(3), g(3)
      REAL(dbl) :: b1(3), b2(3), b3(3)
      LOGICAL :: tqk

! ... end of declarations
!  ----------------------------------------------

! ... mpime = processor number, starting from 0

! ... evaluate cutoffs in reciprocal space and the required mesh size
      sqgc  = sqrt(gcut)
      nr     = int(sqgc) + 2      ! nr   = mesh size parameter

! ... reciprocal lattice generators
      call recips(a1, a2, a3, b1, b2, b3)
      b1 = b1 * alat
      b2 = b2 * alat
      b3 = b3 * alat

! ... verify that, for G<gcut, coordinates never exceed nr
! ... (increase nr if needed)
      CALL vec_prod(c,b1,b2)
      nr3tab=nint(2.0d0*sqgc/abs(dot_prod(c,b3))*vec_mod(c))
      CALL vec_prod(c,b3,b1)
      nr2tab=nint(2.0d0*sqgc/abs(dot_prod(c,b2))*vec_mod(c))
      CALL vec_prod(c,b2,b3)
      nr1tab=nint(2.0d0*sqgc/abs(dot_prod(c,b1))*vec_mod(c))
      nr = max(nr,nr3tab)
      nr = max(nr,nr2tab)
      nr = max(nr,nr1tab)

! ... initialize some variables
      ng     = 0
      nb     = 0

      IF( ALL( qk == 0.0d0 ) ) THEN
        tqk = .FALSE.
      ELSE
        tqk = .TRUE.
      END IF

! *** START OF LOOP ***
! ... calculate moduli of G vectors and the range of indexes where
! ... |G| < gcut 

      DO k = -nr, nr
        IF( MOD( k + nr, nproc ) == mpime ) THEN
          DO j = -nr, nr
            DO i = -nr, nr

              g( 1 ) = REAL(i) * b1(1) + REAL(j) * b2(1) + REAL(k) * b3(1)
              g( 2 ) = REAL(i) * b1(2) + REAL(j) * b2(2) + REAL(k) * b3(2)
              g( 3 ) = REAL(i) * b1(3) + REAL(j) * b2(3) + REAL(k) * b3(3)

! ...         calculate modulus
              IF( tqk ) THEN
                gsq = ( g( 1 ) + qk( 1 ) )**2 + &
                    & ( g( 2 ) + qk( 2 ) )**2 + &
                    & ( g( 3 ) + qk( 3 ) )**2 
              ELSE
                gsq =  g( 1 )**2 + g( 2 )**2 + g( 3 )**2 
              END IF

              IF( gsq < gcut ) THEN
! ...           increase counters
                ng  = ng  + 1
! ...           calculate minimum and maximum indexes
                nb(1) = MAX( nb(1), ABS( i ) )
                nb(2) = MAX( nb(2), ABS( j ) )
                nb(3) = MAX( nb(3), ABS( k ) )
              END IF

            END DO
          END DO
        END IF
      END DO

      CALL mp_sum( ng ,group )

! ... the size of the required (3-dimensional) matrix depends on the
! ... minimum and maximum indices
      CALL mp_max( nb, group )

      nr1 = 2 * nb(1) + 1
      nr2 = 2 * nb(2) + 1
      nr3 = 2 * nb(3) + 1

      RETURN

   CONTAINS

!  ----------------------------------------------
!  ----------------------------------------------
      FUNCTION dot_prod(a,b)

!  this function calculates the dot product of two vectors
!  ----------------------------------------------
      
      REAL(dbl) dot_prod
! ... declare function arguments
      REAL(dbl) a(3),b(3)     

! ... evaluate dot product
      dot_prod=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)

      RETURN
      END FUNCTION dot_prod

!  ----------------------------------------------
!  ----------------------------------------------
      FUNCTION vec_mod(a)

!  this function calculates the norm of a vector
!  ----------------------------------------------

      REAL(dbl) vec_mod
! ... declare function argument
      REAL(dbl) a(3)

! ... evaluate norm
      vec_mod=sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))

      RETURN
      END FUNCTION vec_mod

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE vec_prod(c,a,b)
   
!  this subroutine calculates the vector (cross) product of vectors
!  a,b and stores the result in vector c
!  ----------------------------------------------

! ... declare subroutine arguments
      REAL(dbl) a(3),b(3),c(3)

! ... evaluate cross product
      c(1) = a(2)*b(3)-a(3)*b(2)
      c(2) = a(3)*b(1)-a(1)*b(3)
      c(3) = a(1)*b(2)-a(2)*b(1)

      RETURN 
      END SUBROUTINE vec_prod

   END  SUBROUTINE
