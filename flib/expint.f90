!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
      FUNCTION EXPINT(n, x)
!-----------------------------------------------------------------------
!
! Evaluates the exponential integral E_n(x)
! Parameters: maxit is the maximum allowed number of iterations,
! eps is the desired relative error, not smaller than the machine precision,
! big is a number near the largest representable floating-point number,
! Inspired from Numerical Recipes
! 
      USE kinds, ONLY : DP
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n
      REAL(DP), INTENT(IN) :: x
      REAL(DP) :: expint
      INTEGER, parameter :: maxit=200
      REAL(DP), parameter :: eps=1E-12_DP, big=huge(x)*eps
      REAL(DP), parameter :: euler = 0.577215664901532860606512_DP
!     EPS=1E-9, FPMIN=1E-30

      INTEGER :: i, nm1, k
      REAL(DP) :: a,b,c,d,del,fact,h,iarsum

      IF (.NOT. ((n >= 0).AND.(x >= 0.0).AND.((x > 0.0).OR.(n > 1)))) THEN
         CALL errore('expint','bad arguments', 1)
      END IF

      IF (n == 0) THEN
         expint = exp(-x)/x
         RETURN
      END IF
      nm1 = n-1
      IF (x == 0.0_DP) THEN
         expint = 1.0_DP/nm1
      ELSE IF (x > 1.0_DP) THEN
         b = x+n
         c = big
         d = 1.0_DP/b
         h = d
         DO i=1,maxit
            a = -i*(nm1+i)
            b = b+2.0_DP
            d = 1.0_DP/(a*d+b)
            c = b+a/c
            del = c*d
            h = h*del
            IF (ABS(del-1.0_DP) <= EPS) EXIT
         END DO
         IF (i > maxit) CALL errore('expint','continued fraction failed',1)
         expint = h*EXP(-x)
      ELSE
         IF (nm1 /= 0) THEN
            expint = 1.0_DP/nm1
         ELSE
            expint = -LOG(x)-euler
         END IF
         fact = 1.0_DP
         do i=1,maxit
            fact = -fact*x/i
            IF (i /= nm1) THEN
               del = -fact/(i-nm1)
            ELSE

               iarsum = 0.0_DP
               do k=1,nm1
                  iarsum = iarsum + 1.0_DP/k
               end do

               del = fact*(-LOG(x)-euler+iarsum)
!               del = fact*(-LOG(x)-euler+sum(1.0_DP/arth(1,1,nm1)))
            END IF
            expint = expint+del
            IF (ABS(del) < ABS(expint)*eps) EXIT
         END DO
         IF (i > maxit) CALL errore('expint','series failed',1)
      END IF

      END FUNCTION EXPINT

! -------------------------------------------------------------------
