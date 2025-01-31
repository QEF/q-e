!
! Copyright (C) 2003-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file initiated by Carlo Cavazzoni 2020
!
! Purpose: collect miscellaneous subroutines to help dealing with 
!          accelerator devices

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

SUBROUTINE MYZGERC ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
!     .. Scalar Arguments ..
    COMPLEX*16, INTENT(IN) :: ALPHA
    INTEGER,    INTENT(IN) :: INCX, INCY, LDA, M, N
!     .. Array Arguments ..
    COMPLEX*16 :: A( LDA, * ), X( * ), Y( * )
#if defined(__CUDA)
    attributes(device) :: A, X, Y
    CALL cublasZgerc( M, N, ALPHA, X, INCX, Y, INCY, A, LDA)
#else
    CALL ZGERC  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
#endif

END SUBROUTINE MYZGERC

!=----------------------------------------------------------------------------=!

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

SUBROUTINE MYZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
      COMPLEX*16, INTENT(IN) :: ALPHA,BETA
      INTEGER, INTENT(IN) :: INCX,INCY,LDA,M,N
      CHARACTER*1, INTENT(IN) :: TRANS
      COMPLEX*16 :: A(LDA,*),X(*),Y(*)
#if defined(__CUDA)
    attributes(device) :: A, X, Y
    CALL cublaszgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#else
    CALL zgemv(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#endif
END SUBROUTINE MYZGEMV

!=----------------------------------------------------------------------------=

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
!$acc routine( MYDDOT_VECTOR_GPU ) vector
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine MYDTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) 
#if defined(__CUDA)
USE cublas
#endif
implicit none
  character*1 :: side, uplo, transa, diag 
  integer :: m, n, lda, ldb 
  DOUBLE PRECISION :: alpha 
  DOUBLE PRECISION, dimension(lda, *) :: a 
  DOUBLE PRECISION, dimension(ldb, *) :: b 
#if defined(__CUDA)
  attributes(device) :: a, b 
  call cublasDTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)  
#else
  call DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)  
#endif
  return
end subroutine MYDTRSM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE COMPLEX function MYZDOTC(n, zx, incx, zy, incy) 
#if defined(__CUDA)
USE cublas
#endif
implicit none
  integer :: n, incx, incy 
  DOUBLE COMPLEX, dimension(*) :: zx, zy
#if defined(__CUDA)
  attributes(device) :: zx, zy 
  MYZDOTC = cublasZDOTC(n, zx, incx, zy, incy)  
#else
  DOUBLE COMPLEX, EXTERNAL :: ZDOTC
  MYZDOTC = ZDOTC(n, zx, incx, zy, incy)  
#endif
  return
end function MYZDOTC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE COMPLEX function MYZDOTC_VECTOR_GPU(n, zx, zy) 
#if defined(__CUDA)
!$acc routine(MYZDOTC_VECTOR_GPU) vector
#endif
implicit none
  integer :: n
  DOUBLE COMPLEX, dimension(*) :: zx, zy
#if defined(__CUDA)
  attributes(device) :: zx, zy 
  !civn: from here lapack source code. code for unequal increments 
  !      or equal increments not equal to 1 NOT implemented
  COMPLEX*16 ztemp
  INTEGER i
  INTRINSIC dconjg
  ztemp = (0.0d0,0.0d0)
  MYZDOTC_VECTOR_GPU = (0.0d0,0.0d0)
  IF (n.LE.0) RETURN
  !$acc loop vector reduction(+:ztemp) 
  DO i = 1,n
     ztemp = ztemp + dconjg(zx(i))*zy(i)
  END DO
  MYZDOTC_VECTOR_GPU = ztemp
#else
  DOUBLE COMPLEX, EXTERNAL :: ZDOTC
  MYZDOTC_VECTOR_GPU = ZDOTC(n, zx, 1, zy, 1)  
#endif
  return
end function MYZDOTC_VECTOR_GPU 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZSWAP(n, zx, incx, zy, incy) 
#if defined(__CUDA)
USE cublas
#endif
implicit none
  integer :: n, incx, incy 
  DOUBLE COMPLEX, dimension(*) :: zx, zy
#if defined(__CUDA)
  attributes(device) :: zx, zy 
  CALL cublasZSWAP(n, zx, incx, zy, incy)  
#else
  CALL ZSWAP(n, zx, incx, zy, incy)  
#endif
  return
END SUBROUTINE MYZSWAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZSWAP_VECTOR_GPU(n, zx, zy) 
#if defined(__CUDA)
!$acc routine(MYZSWAP_VECTOR_GPU) vector
#endif
implicit none
  integer :: n
  DOUBLE COMPLEX, dimension(*) :: zx, zy
#if defined(__CUDA)
  attributes(device) :: zx, zy 
  complex*16 ztemp
  integer i
  IF (n.LE.0) RETURN
  !$acc loop vector private(ztemp)
  DO i = 1,n
     ztemp = zx(i)
     zx(i) = zy(i)
     zy(i) = ztemp
  END DO 
#else
  CALL ZSWAP(n, zx, 1, zy, 1)  
#endif
  return
END SUBROUTINE MYZSWAP_VECTOR_GPU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZCOPY(n, zx, incx, zy, incy)
#if defined(__CUDA)
USE cublas
#endif
IMPLICIT NONE 
  INTEGER :: n, incx, incy
  DOUBLE COMPLEX, dimension(*) :: zx, zy
#if defined(__CUDA)
  attributes(device) :: zx, zy 
  CALL cublasZCOPY(n, zx, incx, zy, incy)  
#else
  CALL ZCOPY(n, zx, incx, zy, incy)  
#endif
  RETURN
END SUBROUTINE MYZCOPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZAXPY(n, za, zx, incx, zy, incy)
#if defined(__CUDA)
USE cublas
#endif
IMPLICIT NONE 
  INTEGER :: n, incx, incy
  DOUBLE COMPLEX :: za  
  DOUBLE COMPLEX, dimension(*) :: zx, zy
#if defined(__CUDA)
  attributes(device) :: zx, zy 
  CALL cublasZAXPY(n, za, zx, incx, zy, incy)  
#else
  CALL ZAXPY(n, za, zx, incx, zy, incy)  
#endif
  RETURN
END SUBROUTINE MYZAXPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZDSCAL(n, da, zx, incx)
#if defined(__CUDA)
USE cublas
#endif
IMPLICIT NONE 
  INTEGER :: n, incx
  DOUBLE PRECISION :: da
  DOUBLE COMPLEX, dimension(*) :: zx
#if defined(__CUDA)
  attributes(device) :: zx
  CALL cublasZDSCAL(n, da, zx, incx)
#else
  CALL ZDSCAL(n, da, zx, incx)
#endif
  RETURN
END SUBROUTINE MYZDSCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYZSCAL(n, za, zx, incx)
#if defined(__CUDA)
USE cublas
#endif
IMPLICIT NONE 
  INTEGER :: n, incx
  DOUBLE COMPLEX :: za
  DOUBLE COMPLEX, dimension(*) :: zx
#if defined(__CUDA)
  attributes(device) :: zx
  CALL cublasZSCAL(n, za, zx, incx)
#else
  CALL ZSCAL(n, za, zx, incx)
#endif
  RETURN
END SUBROUTINE MYZSCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYDCOPY(n, x, incx, y, incy)
#if defined(__CUDA)
USE cublas
#endif
IMPLICIT NONE
  INTEGER :: n, incx, incy
  DOUBLE PRECISION, INTENT(IN)   :: x(*)
  DOUBLE PRECISION, INTENT(OUT)  :: y(*)
#if defined(__CUDA)
  attributes(device) :: x, y
  call cublasDCOPY(n, x, incx, y, incy)
#else
  call DCOPY(n, x, incx, y, incy)
#endif
  RETURN
END SUBROUTINE MYDCOPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYDAXPY(n, a, x, incx, y, incy)
#if defined(__CUDA)
USE cublas
#endif
IMPLICIT NONE
  INTEGER :: n, incx, incy
  DOUBLE PRECISION, INTENT(IN)  :: a
  DOUBLE PRECISION, INTENT(IN)  :: x(*) 
  DOUBLE PRECISION, INTENT(OUT) :: y(*) 
#if defined(__CUDA)
  attributes(device) :: x, y
  call cublasDAXPY( n, a, x, incx, y, incy)
#else
  call DAXPY( n, a, x, incx, y, incy)
#endif
  RETURN
END SUBROUTINE MYDAXPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYDSCAL(n, a, x, incx)
#if defined(__CUDA)
USE cublas
#endif
IMPLICIT NONE
  integer :: n, incx
  DOUBLE PRECISION :: a
  DOUBLE PRECISION, dimension(*)  :: x
#if defined(__CUDA)
  attributes(device) :: x
  call cublasDSCAL(n, a, x, incx)
#else
  call DSCAL(n, a, x, incx)
#endif
  RETURN
END SUBROUTINE MYDSCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYDSWAP(n, dx, incx, dy, incy) 
#if defined(__CUDA)
USE cublas
#endif
implicit none
  integer :: n, incx, incy 
  REAL(8), dimension(*) :: dx, dy
#if defined(__CUDA)
  attributes(device) :: dx, dy 
  CALL cublasDSWAP(n, dx, incx, dy, incy)  
#else
  CALL DSWAP(n, dx, incx, dy, incy)  
#endif
  return
END SUBROUTINE MYDSWAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MYDSWAP_VECTOR_GPU(n, dx, dy) 
#if defined(__CUDA)
!$acc routine(MYDSWAP_VECTOR_GPU) vector
#endif
implicit none
  integer :: n
  DOUBLE PRECISION, dimension(*) :: dx, dy
#if defined(__CUDA)
  attributes(device) :: dx, dy 
  double precision dtemp
  integer i
  IF (n.LE.0) RETURN
  !$acc loop vector private(dtemp)
  DO i = 1,n
     dtemp = dx(i)
     dx(i) = dy(i)
     dy(i) = dtemp
  END DO 
#else
  CALL DSWAP(n, dx, 1, dy, 1)  
#endif
  return
END SUBROUTINE MYDSWAP_VECTOR_GPU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
