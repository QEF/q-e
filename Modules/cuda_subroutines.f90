  !----------------------------------------------
  ! ... this file contains a number of subroutines optionally interfaced
  ! ... to cublas
  !----------------------------------------------

SUBROUTINE cudaDGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
    implicit none
    DOUBLE PRECISION :: ALPHA,BETA
    INTEGER :: INCX,INCY,LDA,M,N
    CHARACTER :: TRANS
    DOUBLE PRECISION :: A(LDA,*),X(*),Y(*)
#if defined(__CUDA)
    attributes(device) :: A, X, Y
#endif
    !
    call DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    !
END SUBROUTINE cudaDGEMV

SUBROUTINE cudaDGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
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

END SUBROUTINE cudaDGER
