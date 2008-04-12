!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine pw_dot(sum_over_nodes,n,m,a,lda,b,ldb,c)
!-----------------------------------------------------------------------
  !
  !  calculate m dot products c_i = real( a^*_ij b_ji )
  !  using half G vectors or half PWs
  !
#include "f_defs.h"
  USE kinds, only: DP
  use gvect, only: gstart
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  implicit none
  ! input
  integer :: n, m, lda, ldb
  character(len=1) sum_over_nodes
  complex(DP) :: a(lda,m), b(ldb,m)
  ! output
  real(DP) :: c(m)
  ! local
  integer i
  real(DP), EXTERNAL :: DDOT
  !
  do i= 1,m
     c(i) = 2.d0*DDOT(2*n,a(1,i),1,b(1,i),1)
     if (gstart==2) c(i) = c(i) - DBLE(a(1,i))*DBLE(b(1,i))
  end do
#ifdef __PARA
  if (sum_over_nodes.eq.'y'.or.sum_over_nodes.eq.'Y') call mp_sum( c, intra_pool_comm )
#endif
  return
end subroutine pw_dot
