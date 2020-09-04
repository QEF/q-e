!Generic cublas subroutineU
!
SUBROUTINE gpu_DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
#if defined(__CUDA)
    use cudafor
    use cublas
#endif
    IMPLICIT NONE
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
END SUBROUTINE gpu_DGEMV
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gpu_DCOPY(n, x_d, incx, y_d, incy)
#if defined(__CUDA)
    use cudafor
    USE cublas
#endif
  IMPLICIT NONE
  integer :: n, incx, incy
  DOUBLE PRECISION, INTENT(IN)   :: x_d(*)
  DOUBLE PRECISION, INTENT(OUT)  :: y_d(*)
#if defined(__CUDA)
  attributes(device) :: x_d, y_d
#endif
  !
  call DCOPY(n, x_d, incx, y_d, incy)
  !
END SUBROUTINE gpu_DCOPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gpu_DAXPY(n, a_d, x_d, incx, y_d, incy)
#if defined(__CUDA)
    use cudafor
    USE cublas
#endif
  IMPLICIT NONE
  integer :: n, incx, incy
  DOUBLE PRECISION, INTENT(IN)  :: a_d
  DOUBLE PRECISION, INTENT(IN)  :: x_d(*) 
  DOUBLE PRECISION, INTENT(OUT) :: y_d(*) 
#if defined(__CUDA)
  attributes(device) :: x_d, y_d
#endif
  !
  call cublasDAXPY( n, a_d, x_d, incx, y_d, incy)
  !
END SUBROUTINE gpu_DAXPY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gpu_DSCAL(n, a_d, x_d, incx, y_d, incy)
#if defined(__CUDA)
    use cudafor
    USE cublas
#endif
  IMPLICIT NONE
  integer :: n, incx, incy
  DOUBLE PRECISION :: a_d
  DOUBLE PRECISION, dimension(*)  :: x_d, y_d 
#if defined(__CUDA)
  attributes(device) :: x_d, y_d, a_d
#endif
  !
  call DAXPY(n, a_d, x_d, incx, y_d, incy)
  !
END SUBROUTINE gpu_DSCAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION gpu_DAdd(a_d, b_d)
#if defined(__CUDA)
USE cudafor
#endif
  IMPLICIT NONE
  DOUBLE PRECISION :: a_d, b_d
  DOUBLE PRECISION :: gpu_DAdd, a
#if defined(__CUDA)
  attributes(device) :: a_d, b_d
#endif
  a = a_d
  gpu_DAdd= a + b_d
  return
END FUNCTION gpu_DAdd
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION gpu_DMul(a_d, b_d)
#if defined(__CUDA)
USE cudafor
#endif
  IMPLICIT NONE
  DOUBLE PRECISION :: a_d, b_d
  DOUBLE PRECISION :: gpu_DMul, a
#if defined(__CUDA)
  attributes(device) :: a_d, b_d
#endif
  a = a_d
  gpu_DMul= a * b_d
  return
END FUNCTION gpu_DMul
