!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Fri Oct  8 15:46:10 MDT; 1999
!  ----------------------------------------------
!  BEGIN manual
!
!  SUBROUTINE miller2nxh(i,j,k,n1h,n2h,n3h,nr1,nr2,nr3)
!  SUBROUTINE miller2inxh(i,j,k,in1h,in2h,in3h,nr1,nr2,nr3)
!  SUBROUTINE miller2indx(i,j,k,ind1,ind2,ind3,nr1,nr2,nr3)
!  SUBROUTINE inxh2miller(in1h,in2h,in3h,i,j,k,nr1,nr2,nr3)
!  SUBROUTINE inxh2nxh(in1h,in2h,in3h,n1h,n2h,n3h,nr1,nr2,nr3)
!  SUBROUTINE inxh2indx(in1h,in2h,in3h,ind1,ind2,ind3,nr1,nr2,nr3)
!  ----------------------------------------------
!  these routines perform conversions between Miller indices and matrix
!  indices. We use two sets of matrix indices, that we label
!  n1h,n2h,n3h and in1h,in2h,in3h respectively. Miller indices are
!  labelled i,j,k
!  allowed ranges for Miller indices: -nr1<i<nr1, -nr2<j<nr2, -nr3<k<nr3
!  transformations:
!    n1h = i+1,     i>=0     n2h = j+1,     j>=0     n3h = k+1,     k>=0
!          i+1+nr1, i<0            j+1+nr2, j<0            k+1+nr3, k<0
!  n1h,n2h,n3h range from 1 to nr1,nr2,nr3 respectively. Indices derived
!  from different Miller indices do not overlap as long as abs(i)<nr1/2,
!  abs(j)<nr2/2, abs(k)<nr3/2
!    in1h = i+1,    i>=0     in2h = j+1,    j>=0     in3h = k+1,    k>=0
!           -i+nr1, i<0             -j+nr2, j<0             -k+nr3, k<0
!  in1h,in2h,in3h range from 1 to 2*nr1-1, 2*nr2-1, 2*nr3-1
!  respectively. There is 1-to-1 correspondence with Miller indices
!  n1h,n2h,n3h indices are used for Fourier transforms. The mapping from
!  Miller indices to n1h,n2h,n3h introduces a factor exp(m*2*pi*i)=1 in
!  the Fourier transform
!  ----------------------------------------------
!  END manual
!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE miller2nxh(i,j,k,n1h,n2h,n3h,nr1,nr2,nr3)

!  convert Miller indices to n1h, n2h, n3h ones
!  ----------------------------------------------

      IMPLICIT NONE

!  declare subroutine arguments
      INTEGER, INTENT(IN)  :: i,j,k
      INTEGER, INTENT(OUT) :: n1h,n2h,n3h
      INTEGER, INTENT(IN)  :: nr1,nr2,nr3

!  end of declarations
!  ----------------------------------------------

      IF (i .GE. 0) THEN
        n1h = i+1
      ELSE
        n1h = i+1+nr1
      END IF

      IF (j .GE. 0) THEN
        n2h = j+1
      ELSE
        n2h = j+1+nr2
      END IF

      IF (k .GE. 0) THEN
        n3h = k+1
      ELSE
        n3h = k+1+nr3
      END IF

      RETURN
      END

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE miller2inxh(i,j,k,in1h,in2h,in3h,nr1,nr2,nr3)

!  convert Miller indices to in1h, in2h, in3h ones
!  ----------------------------------------------

      IMPLICIT NONE

!  declare subroutine arguments
      INTEGER, INTENT(IN)  :: i,j,k
      INTEGER, INTENT(OUT) :: in1h,in2h,in3h 
      INTEGER, INTENT(IN)  :: nr1,nr2,nr3

!  end of declarations
!  ----------------------------------------------

      IF (i .GE. 0) THEN
        in1h = i+1
      ELSE
        in1h = -i+nr1
      END IF

      IF (j .GE. 0) THEN
        in2h = j+1
      ELSE
        in2h = -j+nr2
      END IF

      IF (k .GE. 0) THEN
        in3h = k+1
      ELSE
        in3h = -k+nr3
      END IF

      RETURN
      END

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE miller2indx(i,j,k,ind1,ind2,ind3,nr1,nr2,nr3)

!  convert Miller indices to n1h, n2h, n3h ones for the opposite G-point
!  (by convention we call them ind1, ind2. ind3)
!  ----------------------------------------------

      IMPLICIT NONE

!  declare subroutine arguments
      INTEGER, INTENT(IN)  :: i,j,k
      INTEGER, INTENT(OUT) :: ind1,ind2,ind3 
      INTEGER, INTENT(IN)  :: nr1,nr2,nr3

!  end of declarations
!  ----------------------------------------------

      IF (i .LE. 0) THEN
        ind1 = -i+1
      ELSE
        ind1 = -i+1+nr1
      END IF

      IF (j .LE. 0) THEN
        ind2 = -j+1
      ELSE
        ind2 = -j+1+nr2
      END IF

      IF (k .LE. 0) THEN
        ind3 = -k+1
      ELSE
        ind3 = -k+1+nr3
      END IF

      RETURN
      END

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE inxh2miller(in1h,in2h,in3h,i,j,k,nr1,nr2,nr3)

!  convert in1h, in2h, in3h indices to Miller ones
!  ----------------------------------------------

      IMPLICIT NONE

!  declare subroutine arguments
      INTEGER, INTENT(OUT) :: i,j,k
      INTEGER, INTENT(IN)  :: in1h,in2h,in3h
      INTEGER, INTENT(IN)  :: nr1,nr2,nr3

!  end of declarations
!  ----------------------------------------------

      IF (in1h .LE. nr1) THEN
        i = in1h-1
      ELSE
        i = -in1h+nr1
      END IF

      IF (in2h .LE. nr2) THEN
        j = in2h-1
      ELSE
        j = -in2h+nr2
      END IF

      IF (in3h .LE. nr3) THEN
        k = in3h-1
      ELSE
        k = -in3h+nr3
      END IF

      RETURN
      END

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE inxh2nxh(in1h,in2h,in3h,n1h,n2h,n3h,nr1,nr2,nr3)

!  convert in1h, in2h, in3h indices to n1h, n2h, n3h ones
!  ----------------------------------------------

      IMPLICIT NONE

!  declare subroutine arguments
      INTEGER, INTENT(OUT) :: n1h,n2h,n3h
      INTEGER, INTENT(IN)  :: in1h,in2h,in3h
      INTEGER, INTENT(IN)  :: nr1,nr2,nr3

!  end of declarations
!  ----------------------------------------------

      IF (in1h .LE. nr1) THEN
        n1h = in1h
      ELSE
        n1h = -in1h+1+2*nr1
      END IF

      IF (in2h .LE. nr2) THEN
        n2h = in2h
      ELSE
        n2h = -in2h+1+2*nr2
      END IF

      IF (in3h .LE. nr3) THEN
        n3h = in3h
      ELSE
        n3h = -in3h+1+2*nr3
      END IF

      RETURN
      END

!  ----------------------------------------------
!  ----------------------------------------------
      SUBROUTINE inxh2indx(in1h,in2h,in3h,ind1,ind2,ind3,nr1,nr2,nr3)

!  convert in1h, in2h, in3h indices to n1h, n2h, n3h ones for the
!  opposite G-point (which we call ind1, ind2. ind3)
!  ----------------------------------------------

      IMPLICIT NONE

!  declare subroutine arguments
      INTEGER, INTENT(OUT) :: ind1,ind2,ind3
      INTEGER, INTENT(IN)  :: in1h,in2h,in3h
      INTEGER, INTENT(IN)  :: nr1,nr2,nr3

!  end of declarations
!  ----------------------------------------------

      IF (in1h .EQ. 1) THEN
        ind1 = 1
      ELSE IF (in1h .LE. nr1) THEN
        ind1 = 2-in1h+nr1
      ELSE
        ind1 = in1h-nr1+1
      END IF

      IF (in2h .EQ. 1) THEN
        ind2 = 1
      ELSE IF (in2h .LE. nr2) THEN
        ind2 = 2-in2h+nr2
      ELSE
        ind2 = in2h-nr2+1
      END IF

      IF (in3h .EQ. 1) THEN
        ind3 = 1
      ELSE IF (in3h .LE. nr3) THEN
        ind3 = 2-in3h+nr3
      ELSE
        ind3 = in3h-nr3+1
      END IF

      RETURN
      END

!  ----------------------------------------------
!  ----------------------------------------------
      FUNCTION miller2gsq(i,j,k,b1,b2,b3)
    
!  calculate squared modulus of G
!  ----------------------------------------------

      USE kinds
      IMPLICIT NONE
      REAL(dbl) miller2gsq
!  declare subroutine arguments
      INTEGER i,j,k
      REAL(dbl) b1(3), b2(3), b3(3)

!  declare other variables
      REAL(dbl) gsq

!  end of declarations
!  ----------------------------------------------

      gsq =       ( REAL(i)*b1(1) + REAL(j)*b2(1) + REAL(k)*b3(1) ) ** 2
      gsq = gsq + ( REAL(i)*b1(2) + REAL(j)*b2(2) + REAL(k)*b3(2) ) ** 2
      gsq = gsq + ( REAL(i)*b1(3) + REAL(j)*b2(3) + REAL(k)*b3(3) ) ** 2

      miller2gsq = gsq

      RETURN
      END
