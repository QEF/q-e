!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!----------------------------------------------------------------------
subroutine invmat (a, a_inv, n)
  !----------------------------------------------------------------------
  ! computes the inverse "a_inv" of matrix "a", both dimensioned (n,n)
  ! matrix "a" is unchanged on output - LAPACK
  !
  use parameters
  implicit none
  integer :: n
  real(kind=DP) :: a (n, n), a_inv (n, n)

  integer, parameter :: lwork = 100
  ! lwork is the dimension of work and ipiv. Must be >= n
  integer :: info, lda, ipiv (lwork)
  ! info=0: inversion was successful
  ! lda   : leading dimension (the same as n)
  ! ipiv  : work space for pivoting
  real(kind=DP) :: work (lwork)

  if (n.gt.lwork) call error ('invmat', 'work space is too small ', n)
  lda = n
  call DCOPY (n * n, a, 1, a_inv, 1)
  call DGETRF (n, n, a_inv, lda, ipiv, info)
  call error ('invmat', 'error in DGETRF', abs (info) )
  call DGETRI (n, a_inv, lda, ipiv, work, lwork, info)
  call error ('invmat', 'error in DGETRI', abs (info) )

  return

end subroutine invmat

