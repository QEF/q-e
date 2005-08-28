!
! Copyright (C) 2002-2005 FPMD-CPV groups
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
!  Last modified: Thu Oct  7 14:55:48 MDT; 1999
!  ----------------------------------------------
      FUNCTION miller2gsqr(i,j,k,b1,b2,b3)
    
!  calculate squared modulus of G
!  ----------------------------------------------

      USE kinds
      IMPLICIT NONE

      REAL(DP) miller2gsqr
!  declare subroutine arguments
      INTEGER i,j,k
      REAL(DP) b1(3), b2(3), b3(3)

!  declare other variables
      REAL(DP) gsq

!  end of declarations
!  ----------------------------------------------

      gsq =       ( DBLE(i)*b1(1) + DBLE(j)*b2(1) + DBLE(k)*b3(1) ) ** 2
      gsq = gsq + ( DBLE(i)*b1(2) + DBLE(j)*b2(2) + DBLE(k)*b3(2) ) ** 2
      gsq = gsq + ( DBLE(i)*b1(3) + DBLE(j)*b2(3) + DBLE(k)*b3(3) ) ** 2

      miller2gsqr = gsq

      RETURN
      END FUNCTION miller2gsqr

