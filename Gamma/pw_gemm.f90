!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

! ccalbec (nkb, npwx, npw, nbnd, vkb, psi, bec) =>
!    pw_gemm ('Y', nkb, nbnd, npw, vkb, npwx, psi, npwx, bec, nkb)
!
!-----------------------------------------------------------------------
subroutine pw_gemm (sum_over_nodes, na, nb, n, a, lda, b, ldb, c, ldc)
  !-----------------------------------------------------------------------
  !
  !   matrix times matrix with summation index running on G-vectors or PWs
  !   c(ij)=real(a(ik)*b(kj)) using half G vectors or half PWs
  !
#include "machine.h"
  use parameters, only: DP
  use gvect, only: gstart
  implicit none
  ! input
  integer :: na, nb, n, lda, ldb, ldc
  character(len=1) sum_over_nodes
  complex(kind=DP) :: a(lda,na), b(ldb,nb)
  ! output
  real(kind=DP) :: c(ldc,nb)

  if (na.eq.0.or.nb.eq.0) return

  call start_clock ('pw_gemm')
  call DGEMM ('C', 'N', na, nb, 2*n, 2.d0, a, 2*lda, b, 2*ldb, 0.d0, c, ldc)
  if (gstart==2) call DGER ( na, nb, -1.d0, a, 2*lda, b, 2*ldb, c, ldc)
#ifdef __PARA
  if (sum_over_nodes.eq.'y'.or.sum_over_nodes.eq.'Y') call reduce (ldc * nb, c)
#endif
  call stop_clock ('pw_gemm')
  return
end subroutine pw_gemm

