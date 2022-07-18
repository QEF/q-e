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

! this is analogus to MYDDOT, but the result is on device
DOUBLE PRECISION FUNCTION MYDDOT_VECTOR_GPU(N,DX,DY)
#if defined(__CUDA)
!$acc routine(MYDDOT_VECTOR_GPU) vector
#endif
    INTEGER, INTENT(IN) :: N
    DOUBLE PRECISION, INTENT(IN) :: DX(*),DY(*)
    DOUBLE PRECISION :: RES
#if defined (__CUDA)
    INTEGER :: I, M, MP1
    RES = 0.0d0
    MYDDOT_VECTOR_GPU = 0.0d0 
    IF (N.LE.0) RETURN
    ! code for unequal increments or equal increments not equal to 1 NOT implemented 
    M = mod(N,5)
    IF(M.NE.0) THEN
      !$acc loop vector reduction(+:RES)
      DO I = 1, M
        RES = RES + DX(I) * DY(I) 
      END DO 
      IF (N.LT.5) THEN
        MYDDOT_VECTOR_GPU = RES
        RETURN
      END IF
    END IF 
    MP1 = M + 1 
    !$acc loop vector reduction(+:RES)
    DO I = MP1, n, 5
      RES = RES + DX(I)*DY(I) + DX(I+1)*DY(I+1) + DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
    END DO 
    MYDDOT_VECTOR_GPU = RES
#else
    DOUBLE PRECISION DDOT
    MYDDOT_VECTOR_GPU = DDOT(N,DX,1,DY,1)
#endif
END FUNCTION MYDDOT_VECTOR_GPU

function MYDDOTv2 (n, dx, incx, dy, incy)
#if defined(__CUDA)
USE cudafor
USE cublas
#endif
implicit none
  DOUBLE PRECISION :: MYDDOTv2
  integer :: n, incx, incy
  DOUBLE PRECISION, dimension(*)  :: dx, dy
#if defined(__CUDA)
  attributes(device) :: dx, dy
  attributes(device) :: MYDDOTv2
  type(cublashandle) :: h
  integer :: ierr
  h = cublasGetHandle()
  ierr = cublasDDOT_v2(h, n, dx, incx, dy, incy, MYDDOTv2)
#else
  DOUBLE PRECISION DDOT
  MYDDOTv2=DDOT(n, dx, incx, dy, incy)
#endif

  return
end function MYDDOTv2

subroutine MYDDOTv3 (n, dx, incx, dy, incy, result)
#if defined(__CUDA)
USE cudafor
USE cublas
#endif
  implicit none
  integer :: n, incx, incy
  DOUBLE PRECISION, dimension(*)  :: dx, dy
  DOUBLE PRECISION :: result
#if defined(__CUDA)
  attributes(device) :: dx, dy
  attributes(device) :: result
  type(cublashandle) :: h
  integer :: ierr
  h = cublasGetHandle()
  ierr = cublasDDOT_v2(h, n, dx, incx, dy, incy, result)
#else
  DOUBLE PRECISION DDOT
  result=DDOT(n, dx, incx, dy, incy)
#endif

  return
end subroutine MYDDOTv3

