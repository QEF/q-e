!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE pw_dot(sum_over_nodes,n,m,a,lda,b,ldb,c)
!-----------------------------------------------------------------------
  !
  !  calculate m dot products c_i = real( a^*_ij b_ji )
  !  using half G vectors or half PWs
  !
  USE kinds, ONLY: DP
  USE gvect, ONLY: gstart
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  IMPLICIT NONE
  ! input
  INTEGER :: n, m, lda, ldb
  CHARACTER(len=1) sum_over_nodes
  COMPLEX(DP) :: a(lda,m), b(ldb,m)
  ! output
  real(DP) :: c(m)
  ! local
  INTEGER i
  real(DP), EXTERNAL :: ddot
  !
  DO i= 1,m
     c(i) = 2.d0*ddot(2*n,a(1,i),1,b(1,i),1)
     IF (gstart==2) c(i) = c(i) - dble(a(1,i))*dble(b(1,i))
  ENDDO
#if defined(__MPI)
  IF (sum_over_nodes=='y'.or.sum_over_nodes=='Y') CALL mp_sum( c, intra_pool_comm )
#endif
  RETURN
END SUBROUTINE pw_dot
