!
! Copyright (C) 2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE matrix_inversion
   !! Contains LAPACK based routines for matrix inversion.
   !
   IMPLICIT NONE
   PRIVATE
   PUBLIC :: invmat

   INTERFACE invmat
      MODULE PROCEDURE invmat_r, invmat_c
   END INTERFACE

   CONTAINS

  SUBROUTINE invmat_r (n, a, a_inv, da)
  !-----------------------------------------------------------------------
  !! Computes the inverse of a n*n real matrix "a" using LAPACK routines.
  !
  !! * if "a_inv" is not present, "a" contains the inverse on output;
  !! * if "a_inv" is present, it contains the inverse on output, "a" is unchanged;
  !! * if "da" is specified and if the matrix is dimensioned 3x3;
  !! * it also returns the determinant in "da".
  !
#if defined(_OPENMP)
  USE omp_lib
#endif
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  REAL(DP), DIMENSION (n,n), INTENT(inout)  :: a
  REAL(DP), DIMENSION (n,n), INTENT(out), OPTIONAL :: a_inv
  REAL(DP), OPTIONAL, INTENT(out) :: da
  !
  INTEGER :: info, lda, lwork
  ! info=0: inversion was successful
  ! lda   : leading dimension (the same as n)
  INTEGER, ALLOCATABLE :: ipiv (:)
  ! ipiv  : work space for pivoting
  REAL(DP), ALLOCATABLE :: work (:)
  ! more work space
  INTEGER, SAVE :: lworkfact = 64
  !
  ! WORKAROUND STARTS ==================================================
  !
  ! Comment taken from v6.1 (https://github.com/fspiga/qe-gpu)
  !
  ! There is a bug in several versions of MKL that will cause an hang in the multithreaded  DGEMM for AVX2.
  ! To avoid the bug, we have two options, set the number of MKL threads to one or force to use AVX instead of AVX2.
  ! To force the single threads, we need to read the current number of threads with   numt=mkl_get_max_threads(), set it
  ! temporarely to one with "call mkl_set_num_threads(1)" and then  resetting it to the original numt at the end of the function.
  ! To force AVX, we can call mkl_cbwr_set(MKL_CBWR_AVX).
  !
  ! There is currently no way to check if MKL is used for LA.
  ! Since the size of the matrix to be inverted in PW is generally small,
  ! we disable OpenMP in the calls below for the time being.
  !
#ifdef _OPENMP
  INTEGER :: num_threads
  num_threads=omp_get_max_threads()
  CALL omp_set_num_threads(1) 
#endif
  ! WORKAROUND ENDS ====================================================
  !
  IF ( PRESENT(da) ) THEN
     IF ( n == 3 ) THEN
        da = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
             a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
             a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
        IF (abs(da) < 1.d-10) CALL errore(' invmat ',' singular matrix ', 1)
     ELSE
        da = 0.0_dp
     ENDIF
  ENDIF
  !
  lda = n
  lwork=64*n
  ALLOCATE(ipiv(n), work(lwork) )
  !
  IF ( PRESENT(a_inv) ) THEN
     a_inv(:,:) = a(:,:)
     CALL dgetrf (n, n, a_inv, lda, ipiv, info)
  ELSE
     CALL dgetrf (n, n, a, lda, ipiv, info)
  END IF
  CALL errore ('invmat', 'error in DGETRF', abs (info) )
  IF ( PRESENT(a_inv) ) THEN
     CALL dgetri (n, a_inv, lda, ipiv, work, lwork, info)
  ELSE
     CALL dgetri (n, a, lda, ipiv, work, lwork, info)
  END IF 
  CALL errore ('invmat', 'error in DGETRI', abs (info) )
  !
  lworkfact = INT (work(1)/n)
  DEALLOCATE ( work, ipiv )
  ! WORKAROUND STARTS ==================================================
  !
  ! ... and now restrore the previous value.
  !
#ifdef _OPENMP
  CALL omp_set_num_threads(num_threads)
#endif
  ! WORKAROUND ENDS ====================================================
  END SUBROUTINE invmat_r

  SUBROUTINE invmat_c (n, a, a_inv, da)
  !-----------------------------------------------------------------------
  !! Computes the inverse of a n*n complex matrix "a" using LAPACK routines.
  !! See \(\texttt{invmat_r}\) for more infos.
  !
#if defined(_OPENMP)
  USE omp_lib
#endif
  USE kinds, ONLY : DP
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  COMPLEX (DP), DIMENSION (n,n), INTENT(INOUT)  :: a
  COMPLEX (DP), OPTIONAL, DIMENSION (n,n), INTENT(OUT) :: a_inv
  COMPLEX (DP), OPTIONAL, INTENT(OUT) :: da
  !
  INTEGER :: info, lda, lwork
  ! info=0: inversion was successful
  ! lda   : leading dimension (the same as n)
  INTEGER, ALLOCATABLE :: ipiv (:)
  ! ipiv  : work space for pivoting
  COMPLEX(DP), ALLOCATABLE :: work (:)
  ! more work space
  INTEGER, SAVE :: lworkfact = 64
  !
  ! WORKAROUND STARTS ==================================================
  !
  ! Comment taken from v6.1 (https://github.com/fspiga/qe-gpu)
  !
  ! There is a bug in several versions of MKL that will cause an hang in the multithreaded  DGEMM for AVX2.
  ! To avoid the bug, we have two options, set the number of MKL threads to one or force to use AVX instead of AVX2.
  ! To force the single threads, we need to read the current number of threads with   numt=mkl_get_max_threads(), set it
  ! temporarely to one with "call mkl_set_num_threads(1)" and then  resetting it to the original numt at the end of the function.
  ! To force AVX, we can call mkl_cbwr_set(MKL_CBWR_AVX).
  !
  ! There is currently no way to check if MKL is used for LA.
  ! Since the size of the matrix to be inverted in PW is generally small,
  ! we disable OpenMP in the calls below for the time being.
  !
#ifdef _OPENMP
  INTEGER :: num_threads
  num_threads=omp_get_max_threads()
  CALL omp_set_num_threads(1) 
#endif
  ! WORKAROUND ENDS ====================================================
  !
  IF ( PRESENT(da) ) THEN
     IF (n == 3) THEN
        da = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
             a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
             a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
        IF (abs(da) < 1.d-10) CALL errore(' invmat ',' singular matrix ', 1)
     ELSE
        da = (0.d0,0.d0)
     ENDIF
  ENDIF
  !
  lda = n
  lwork=64*n
  ALLOCATE(ipiv(n), work(lwork) )
  !
  IF ( PRESENT(a_inv) ) THEN
     a_inv(:,:) = a(:,:)
     CALL zgetrf (n, n, a_inv, lda, ipiv, info)
  ELSE
     CALL zgetrf (n, n, a, lda, ipiv, info)
  END IF
  CALL errore ('invmat', 'error in ZGETRF', abs (info) )
  IF ( PRESENT(a_inv) ) THEN
     CALL zgetri (n, a_inv, lda, ipiv, work, lwork, info)
  ELSE
     CALL zgetri (n, a, lda, ipiv, work, lwork, info)
  END IF
  CALL errore ('invmat', 'error in ZGETRI', abs (info) )
  !
  lworkfact = INT (work(1)/n)
  DEALLOCATE ( work, ipiv )
  !
  ! WORKAROUND STARTS ==================================================
  !
  ! ... and now restrore the previous value.
  !
#ifdef _OPENMP
  CALL omp_set_num_threads(num_threads)
#endif
  ! WORKAROUND ENDS ====================================================
  !
  END SUBROUTINE invmat_c

END MODULE matrix_inversion
