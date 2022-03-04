!
! Copyright (C) 2003-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file initiated by Carlo Cavazzoni 2020
!
! Purpose: collect miscellaneus subroutines to help dealing with 
!          accelerator devices

! In principle this can go away .......
SUBROUTINE MYDGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
!     .. Scalar Arguments ..
    DOUBLE PRECISION ::  ALPHA
    INTEGER          ::   INCX, INCY, LDA, M, N
!     .. Array Arguments ..
    DOUBLE PRECISION :: A( LDA, * ), X( * ), Y( * )
#if defined(__CUDA)
    attributes(device) :: A, X, Y
#endif
    CALL DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )

END SUBROUTINE MYDGER

!=----------------------------------------------------------------------------=!

! In principle this can go away .......
SUBROUTINE MYDGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
    CHARACTER*1, INTENT(IN) ::        TRANSA, TRANSB
    INTEGER, INTENT(IN) ::            M, N, K, LDA, LDB, LDC
    DOUBLE PRECISION, INTENT(IN) ::   ALPHA, BETA
    DOUBLE PRECISION  :: A( LDA, * ), B( LDB, * ), C( LDC, * )
#if defined(__CUDA)
    attributes(device) :: A, B, C
    CALL cublasdgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#else
    CALL dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#endif

END SUBROUTINE MYDGEMM

SUBROUTINE MYZGEMM( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
    CHARACTER*1, INTENT(IN) ::        TRANSA, TRANSB
    INTEGER, INTENT(IN) ::            M, N, K, LDA, LDB, LDC
    COMPLEX*16, INTENT(IN) ::   ALPHA, BETA
    COMPLEX*16  :: A( LDA, * ), B( LDB, * ), C( LDC, * )
#if defined(__CUDA)
    attributes(device) :: A, B, C
    CALL cublaszgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#else
    CALL zgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
#endif

END SUBROUTINE MYZGEMM

! In principle this can go away .......
SUBROUTINE MYDGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
      DOUBLE PRECISION, INTENT(IN) :: ALPHA,BETA
      INTEGER, INTENT(IN) :: INCX,INCY,LDA,M,N
      CHARACTER*1, INTENT(IN) :: TRANS
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
#if defined(__CUDA)
    attributes(device) :: A, X, Y
    CALL cublasdgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#else
    CALL dgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#endif
END SUBROUTINE MYDGEMV

!=----------------------------------------------------------------------------=

! In principle this can go away .......
DOUBLE PRECISION FUNCTION MYDDOT(N,DX,INCX,DY,INCY)
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
    INTEGER, INTENT(IN) :: INCX,INCY,N
    DOUBLE PRECISION, INTENT(IN) :: DX(*),DY(*)
#if defined(__CUDA)
    attributes(device) :: DX, DY
    MYDDOT = CUBLASDDOT(N,DX,INCX,DY,INCY)
#else
    DOUBLE PRECISION DDOT
    MYDDOT = DDOT(N,DX,INCX,DY,INCY)
#endif

END FUNCTION MYDDOT

