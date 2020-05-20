!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE upf_invmat

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: invmat

CONTAINS

  SUBROUTINE invmat (n, a, a_inv)
  !-----------------------------------------------------------------------
  ! computes the inverse of a n*n real matrix "a" using LAPACK routines 
  ! "a_inv" contains the inverse on output; "a" is unchanged
  !
  USE upf_kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  REAL(DP), DIMENSION (n,n), INTENT(inout)  :: a
  REAL(DP), DIMENSION (n,n), INTENT(out)    :: a_inv
  !
  INTEGER :: info, lda, lwork
  ! info=0: inversion was successful
  ! lda   : leading dimension (the same as n)
  INTEGER, ALLOCATABLE :: ipiv (:)
  ! ipiv  : work space for pivoting
  REAL(DP), ALLOCATABLE :: work (:)
  ! more work space
  ! INTEGER, SAVE :: lworkfact = 64
  !
  lda = n
  lwork=64*n
  ALLOCATE(ipiv(n), work(lwork) )
  !
  a_inv(:,:) = a(:,:)
  CALL dgetrf (n, n, a_inv, lda, ipiv, info)
  CALL upf_error ('invmat', 'error in DGETRF', abs (info) )
  CALL dgetri (n, a_inv, lda, ipiv, work, lwork, info)
  CALL upf_error ('invmat', 'error in DGETRI', abs (info) )
  !
  ! lworkfact = INT (work(1)/n)
  DEALLOCATE ( work, ipiv )

  END SUBROUTINE invmat

END MODULE upf_invmat

