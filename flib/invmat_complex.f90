!
! Copyright (C) 2004 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE invmat_complex (n, a, a_inv, da)
  !-----------------------------------------------------------------------
  ! computes the inverse "a_inv" of a complex matrix "a", both 
  ! dimensioned (n,n). If the matrix is dimensioned 3x3, it also computes 
  ! determinant "da". Matrix "a" is unchanged on output - LAPACK
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER :: n
  COMPLEX (DP), DIMENSION (n,n) :: a, a_inv
  COMPLEX (DP) :: da
  !
  INTEGER :: info, lda, lwork, ipiv (n)
  ! info=0: inversion was successful
  ! lda   : leading dimension (the same as n)
  ! ipiv  : work space for pivoting (assumed of length lwork=n)
  COMPLEX (DP) :: work (n) 
  ! more work space
  !
  lda = n
  lwork=n
  !
  a_inv(:,:) = a(:,:)
  !
  CALL zgetrf (n, n, a_inv, lda, ipiv, info)
  CALL errore ('invmat', 'error in ZGETRF', abs (info) )
  CALL zgetri (n, a_inv, lda, ipiv, work, lwork, info)
  CALL errore ('invmat', 'error in ZGETRI', abs (info) )
  !
  IF (n == 3) THEN
     da = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
          a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
          a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
     IF (ABS(da) < 1.d-10) CALL errore(' invmat ',' singular matrix ', 1)
  ELSE
     da = (0.d0,0.d0)
  END IF
  !
  RETURN
  !
END SUBROUTINE invmat_complex
