!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE zvscal(n,lda,m,v,zin,zout)
  IMPLICIT NONE
  INTEGER :: n, lda, m
  real(8) :: v(n), zin(2,lda,m), zout(2,lda,m)
  INTEGER :: i,j
  !
  DO j = 1,m
     DO i = 1,n
        zout(1,i,j) = zin(1,i,j)*v(i)
        zout(2,i,j) = zin(2,i,j)*v(i)
     ENDDO
  ENDDO
  !
  RETURN
END SUBROUTINE zvscal
