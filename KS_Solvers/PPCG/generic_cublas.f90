subroutine gpu_DGEMM (transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
USE cublas
implicit none
  character*1 transa, transb
  integer :: m, n, k, lda, ldb, ldc
  DOUBLE PRECISION :: alpha, beta 
  DOUBLE PRECISION, device, dimension(lda, *)  :: A
  DOUBLE PRECISION, device, dimension(ldb, *)  :: B
  DOUBLE PRECISION, device, dimension(ldc, *)  :: C
  call cublasDGEMM(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  return
end subroutine gpu_DGEMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gpu_ZGEMM (transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
USE cublas
implicit none
  character*1 transa, transb
  integer :: m, n, k, lda, ldb, ldc
  DOUBLE COMPLEX :: alpha, beta 
  DOUBLE COMPLEX, device, dimension(lda, *)  :: A
  DOUBLE COMPLEX, device, dimension(ldb, *)  :: B
  DOUBLE COMPLEX, device, dimension(ldc, *)  :: C
  call cublasZGEMM(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc)
  return
end subroutine gpu_ZGEMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gpu_DGER (m, n, alpha, x, incx, y, incy, a, lda)
USE cublas
implicit none
  integer :: m, n, lda, incx, incy
  DOUBLE PRECISION :: alpha
  DOUBLE PRECISION, device, dimension(lda, *)  :: A
  DOUBLE PRECISION, device, dimension(*)  :: x, y 
  call cublasDGER(m, n, alpha, x, incx, y, incy, a, lda)  
  return
end subroutine gpu_DGER 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function gpu_DDOT (n, dx, incx, dy, incy)
USE cublas
implicit none
  DOUBLE PRECISION :: gpu_DDOT
  integer :: n, incx, incy
  DOUBLE PRECISION, device, dimension(*)  :: dx, dy 
  gpu_DDOT=cublasDDOT(n, dx, incx, dy, incy)  
  return
end function gpu_DDOT 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gpu_DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) 
USE cublas
implicit none
  character*1 :: side, uplo, transa, diag 
  integer :: m, n, lda, ldb 
  DOUBLE PRECISION :: alpha 
  DOUBLE PRECISION, device, dimension(lda, *) :: a 
  DOUBLE PRECISION, device, dimension(ldb, *) :: b 
  call cublasDTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)  
  return
end subroutine gpu_DTRSM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gpu_threaded_memset(array, val, length)
  !
  USE cudafor
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: length
  REAL(DP), device, INTENT(OUT) :: array(length)
  REAL(DP), INTENT(IN) :: val
  !
  INTEGER :: i
  !
  IF (length<=0) RETURN
  !
!$cuf kernel do(1)  
  DO i=1, length
     array(i) = val
  ENDDO
  !
END SUBROUTINE gpu_threaded_memset
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gpu_threaded_assign(array_out, array_in, kdimx, nact, use_idx, idx, bgrp_root_only)
  !  assign (copy) a complex array in a threaded way
  !  array_out( 1:kdimx, 1:nact ) = array_in( 1:kdimx, 1:nact )       or
  !  array_out( 1:kdimx, 1:nact ) = array_in( 1:kdimx, idx(1:nact) )
  !  if the index array idx is given
  !  if  bgrp_root_only is present and .true. the assignement is made only by the 
  !  MPI root process of the bgrp and array_out is zeroed otherwise
  USE cudafor
  USE util_param,   ONLY : DP
  USE mp_bands_util,      ONLY : root_bgrp_id, nbgrp, my_bgrp_id
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: kdimx, nact
  COMPLEX(DP), device, INTENT(OUT) :: array_out( kdimx, nact )
  COMPLEX(DP), device, INTENT(IN)  :: array_in ( kdimx, * )
  INTEGER, device, INTENT(IN) :: idx( * )
  LOGICAL, INTENT(IN) :: bgrp_root_only
  LOGICAL, INTENT(IN) :: use_idx
  !
  INTEGER, PARAMETER :: blocksz = 256
  INTEGER :: nblock
  
  INTEGER :: i, j
  !
  IF (kdimx <=0 .OR. nact<= 0) RETURN
  !
  IF (bgrp_root_only .AND. ( my_bgrp_id /= root_bgrp_id ) ) THEN
     call threaded_memset( array_out, 0.d0, 2*kdimx*nact )
     RETURN
  END IF
  
  nblock = (kdimx - 1)/blocksz  + 1
  
  IF (use_idx ) THEN
!$cuf kernel do(2)
     DO i=1, nact 
       DO j=1,nblock
        array_out(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) = array_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), idx( i ) ) 
       ENDDO 
     ENDDO
  ELSE
!$cuf kernel do(2)
     DO i=1, nact 
       DO j=1,nblock
        array_out(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) = array_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) 
       ENDDO 
     ENDDO
  END IF
  !
END SUBROUTINE gpu_threaded_assign
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE gpu_threaded_backassign(array_out, idx, array_in, kdimx, nact, use_a2, a2_in )
  !  assign (copy) a complex array in a threaded way
  !  array_out( 1:kdimx, idx(1:nact) ) = array_in( 1:kdimx, 1:nact )      or
  !  array_out( 1:kdimx, idx(1:nact) ) = array_in( 1:kdimx, 1:nact )  + a2_in( 1:kdimx, idx(1:nact) (
  !  if a2_in is present
  !  the index array idx is mandatory otherwise one could use previous routine)
  USE cudafor
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)      :: kdimx, nact
  COMPLEX(DP), device, INTENT(INOUT) :: array_out( kdimx, * ) ! we don't want to mess with un referenced columns
  COMPLEX(DP), device, INTENT(IN)  :: array_in ( kdimx, nact )
  COMPLEX(DP), device, INTENT(IN)  :: a2_in ( kdimx, * )
  INTEGER, device, INTENT(IN)      :: idx( * )
  LOGICAL, INTENT(IN) :: use_a2
  !
  INTEGER, PARAMETER :: blocksz = 256
  INTEGER :: nblock

  INTEGER :: i, j
  !
  IF (kdimx <=0 .OR. nact<= 0) RETURN
  !

  nblock = (kdimx - 1)/blocksz  + 1

  IF ( use_a2) THEN
!$cuf kernel do(2)
     DO i=1, nact 
       DO j=1,nblock
        array_out(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), idx( i ) ) = &
            array_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i )     + & 
                a2_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), idx( i ) ) 
       ENDDO 
     ENDDO
  ELSE
!$cuf kernel do(2)
     DO i=1, nact 
       DO j=1,nblock
        array_out(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), idx( i ) ) = array_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) 
       ENDDO 
     ENDDO
  END IF
  !
END SUBROUTINE gpu_threaded_backassign
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
