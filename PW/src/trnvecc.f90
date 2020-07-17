!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE trnvecc( u, at, bg, iflg )
  !-----------------------------------------------------------------------
  !! Transforms a COMPLEX vector in real space (like a displacement)
  !! from crystal to cartesian axis (iflag.gt.0) and viceversa (iflag.le.0
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER :: iflg
  !! input: gives the versus of the transformatio
  REAL(DP) :: at(3,3)
  !! input: the direct lattice vectors
  REAL(DP) :: bg(3,3)
  !! input: the reciprocal lattice vectors
  COMPLEX(DP) :: u(3)
  !! I/O: the vector to transform
  !
  ! ... local variables
  !
  INTEGER :: i, k
  ! counter on polarizations
  COMPLEX(DP) :: wrk(3)
  ! auxiliary variable
  !
  !
  IF (iflg > 0) THEN
     !
     ! forward transformation :
     !
     DO i = 1, 3
        wrk(i) = u(i)
     ENDDO
     !
     DO i = 1, 3
        u(i) = 0.d0
        DO k = 1, 3
           u(i) = u(i) + wrk(k) * at(i,k)
        ENDDO
     ENDDO
  ELSE
     !
     ! backward transformation :
     !
     DO i = 1, 3
        wrk(i) = 0.d0
        DO k = 1, 3
           wrk(i) = wrk(i) + u(k) * bg(k, i)
        ENDDO
     ENDDO
     !
     DO i = 1, 3
        u(i) = wrk(i)
     ENDDO
     !
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE trnvecc
