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


!=----------------------------------------------------------------------------=!

MODULE device_helper
   USE util_param, ONLY : DP
#if defined(__CUDA)
   USE cudafor
#endif
   IMPLICIT NONE
   SAVE
   PRIVATE

   INTERFACE sync_to_device
      MODULE PROCEDURE sync_to_device_c2d, sync_to_device_r2d
   END INTERFACE
   INTERFACE sync_to_host
      MODULE PROCEDURE sync_to_host_c2d, sync_to_host_r2d
   END INTERFACE

   PUBLIC :: sync_to_device
   PUBLIC :: sync_to_host

CONTAINS

   SUBROUTINE sync_to_device_c2d( h, d )
      COMPLEX(DP), INTENT(IN) :: d(:,:)
      COMPLEX(DP), INTENT(OUT) :: h(:,:)
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) :: d
      h = d
#endif
   END SUBROUTINE
   SUBROUTINE sync_to_host_c2d( h, d )
      COMPLEX(DP), INTENT(OUT) :: d(:,:)
      COMPLEX(DP), INTENT(IN) :: h(:,:)
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) :: d
      d = h
#endif
   END SUBROUTINE
   SUBROUTINE sync_to_device_r2d( h, d )
      REAL(DP), INTENT(IN) :: d(:,:)
      REAL(DP), INTENT(OUT) :: h(:,:)
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) :: d
      h = d
#endif
   END SUBROUTINE
   SUBROUTINE sync_to_host_r2d( h, d )
      REAL(DP), INTENT(OUT) :: d(:,:)
      REAL(DP), INTENT(IN) :: h(:,:)
#if defined(__CUDA)
      ATTRIBUTES(DEVICE) :: d
      d = h
#endif
   END SUBROUTINE

END MODULE

!=----------------------------------------------------------------------------=!

SUBROUTINE qe_device_sync()
#if defined(__CUDA)
   USE cudafor
#endif
   INTEGER :: info
#if defined (__CUDA)
   info = cudaDeviceSynchronize()
   IF( info /= 0 ) CALL errore('qe_sync',' error ',ABS(info))
#endif
   RETURN
END SUBROUTINE

!=----------------------------------------------------------------------------=!

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
