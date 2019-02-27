!
! Copyright (C) 2002-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Utility functions to perform memcpy and memset on the device with CUDA Fortran
! cuf_memXXX contains a CUF KERNEL to perform the selected operation
! cu_memsync are wrappers for cuda_memcpy functions
!
MODULE cuda_util
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  PUBLIC :: cuf_memcpy, cuf_memset, cu_memsync
  !
  INTERFACE cuf_memcpy
    MODULE PROCEDURE &
      cuf_memcpy_r1d, &
      cuf_memcpy_r2d, &
      cuf_memcpy_r3d, &
      cuf_memcpy_c1d, &
      cuf_memcpy_c2d, &
      cuf_memcpy_c3d
  END INTERFACE
  !
  INTERFACE cuf_memset
    MODULE PROCEDURE &
      cuf_memset_r1d, &
      cuf_memset_r2d, &
      cuf_memset_r3d, &
      cuf_memset_c1d, &
      cuf_memset_c2d, &
      cuf_memset_c3d
  END INTERFACE
  !
  INTERFACE cu_memsync
    MODULE PROCEDURE &
      h2d_memsync_r1d, &
      h2d_memsync_r2d, &
      h2d_memsync_r3d, &
      h2d_memsync_c1d, &
      h2d_memsync_c2d, &
      h2d_memsync_c3d
#if defined(__CUDA)
    MODULE PROCEDURE &
      d2h_memsync_r1d, &
      d2h_memsync_r2d, &
      d2h_memsync_r3d, &
      d2h_memsync_c1d, &
      d2h_memsync_c2d, &
      d2h_memsync_c3d
#endif
  END INTERFACE
  !
  CONTAINS
  !
  SUBROUTINE cuf_memcpy_r1d(array_out, array_in, range1 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:)
    REAL(DP), INTENT(IN) :: array_in(:)
    INTEGER, INTENT(IN) ::  range1(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    !
    d1s = range1(1)
    d1e = range1(2)
    !
    !$cuf kernel do(1)
    DO i1 = d1s, d1e
       array_out(i1 ) = array_in(i1 )
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_r1d
  !
  SUBROUTINE cuf_memcpy_r2d(array_out, array_in, range1, range2 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:)
    REAL(DP), INTENT(IN) :: array_in(:,:)
    INTEGER, INTENT(IN) ::  range1(2), range2(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    !
    !$cuf kernel do(2)
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2 ) = array_in(i1,i2 )
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_r2d
  !
  SUBROUTINE cuf_memcpy_r3d(array_out, array_in, range1, range2, range3 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:,:)
    REAL(DP), INTENT(IN) :: array_in(:,:,:)
    INTEGER, INTENT(IN) ::  range1(2), range2(2), range3(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    INTEGER :: i3, d3s, d3e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    d3s = range3(1)
    d3e = range3(2)
    !
    !$cuf kernel do(3)
    DO i3 = d3s, d3e
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_r3d
  !
  SUBROUTINE cuf_memcpy_c1d(array_out, array_in, range1 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:)
    COMPLEX(DP), INTENT(IN) :: array_in(:)
    INTEGER, INTENT(IN) ::  range1(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    !
    d1s = range1(1)
    d1e = range1(2)
    !
    !$cuf kernel do(1)
    DO i1 = d1s, d1e
       array_out(i1 ) = array_in(i1 )
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_c1d
  !
  SUBROUTINE cuf_memcpy_c2d(array_out, array_in, range1, range2 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:)
    COMPLEX(DP), INTENT(IN) :: array_in(:,:)
    INTEGER, INTENT(IN) ::  range1(2), range2(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    !
    !$cuf kernel do(2)
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2 ) = array_in(i1,i2 )
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_c2d
  !
  SUBROUTINE cuf_memcpy_c3d(array_out, array_in, range1, range2, range3 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:,:)
    COMPLEX(DP), INTENT(IN) :: array_in(:,:,:)
    INTEGER, INTENT(IN) ::  range1(2), range2(2), range3(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out, array_in
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    INTEGER :: i3, d3s, d3e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    d3s = range3(1)
    d3e = range3(2)
    !
    !$cuf kernel do(3)
    DO i3 = d3s, d3e
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memcpy_c3d
  !
  !
  SUBROUTINE cuf_memset_r1d(array_out, val, range1 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:)
    REAL(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    !
    d1s = range1(1)
    d1e = range1(2)
    !
    !$cuf kernel do(1)
    DO i1 = d1s, d1e
       array_out(i1 ) = val
    ENDDO
    !
  END SUBROUTINE cuf_memset_r1d
  !
  SUBROUTINE cuf_memset_r2d(array_out, val, range1, range2 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:)
    REAL(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2), range2(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    !
    !$cuf kernel do(2)
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2 ) = val
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memset_r2d
  !
  SUBROUTINE cuf_memset_r3d(array_out, val, range1, range2, range3 )
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:,:)
    REAL(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2), range2(2), range3(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    INTEGER :: i3, d3s, d3e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    d3s = range3(1)
    d3e = range3(2)
    !
    !$cuf kernel do(3)
    DO i3 = d3s, d3e
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2,i3 ) = val
    ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memset_r3d
  !
  SUBROUTINE cuf_memset_c1d(array_out, val, range1 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:)
    COMPLEX(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    !
    d1s = range1(1)
    d1e = range1(2)
    !
    !$cuf kernel do(1)
    DO i1 = d1s, d1e
       array_out(i1 ) = val
    ENDDO
    !
  END SUBROUTINE cuf_memset_c1d
  !
  SUBROUTINE cuf_memset_c2d(array_out, val, range1, range2 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:)
    COMPLEX(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2), range2(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    !
    !$cuf kernel do(2)
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2 ) = val
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memset_c2d
  !
  SUBROUTINE cuf_memset_c3d(array_out, val, range1, range2, range3 )
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:,:)
    COMPLEX(DP), INTENT(IN) :: val
    INTEGER, INTENT(IN) ::  range1(2), range2(2), range3(2)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: i1, d1s, d1e
    INTEGER :: i2, d2s, d2e
    INTEGER :: i3, d3s, d3e
    !
    d1s = range1(1)
    d1e = range1(2)
    d2s = range2(1)
    d2e = range2(2)
    d3s = range3(1)
    d3e = range3(2)
    !
    !$cuf kernel do(3)
    DO i3 = d3s, d3e
    DO i2 = d2s, d2e
    DO i1 = d1s, d1e
       array_out(i1,i2,i3 ) = val
    ENDDO
    ENDDO
    ENDDO
    !
  END SUBROUTINE cuf_memset_c3d
  !
  !
  SUBROUTINE h2d_memsync_r1d(array_out, array_in, range1 )
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:)
    REAL(DP), INTENT(IN)    :: array_in(:)
    INTEGER, INTENT(IN) ::  range1(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyHostToDevice )
#else
    array_out(d1_start:d1_start+d1_size-1) = array_in(d1_start:d1_start+d1_size-1)
#endif
    !
  END SUBROUTINE h2d_memsync_r1d
  !
  SUBROUTINE h2d_memsync_r2d(array_out, array_in, range1, range2 )
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:)
    REAL(DP), INTENT(IN)    :: array_in(:,:)
    INTEGER, INTENT(IN) ::  range1(3), range2(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: d2_start, d2_size, d2_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    d2_start = range2(1)
    d2_size = range2(2)
    d2_ld = range2(3)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else

    !pinned_buffer(1:nbase, n_start:n_end) = sc_d( 1:nbase, n_start:n_end )
    !ierr = cudaMemcpy2D( pinned_buffer(1, n_start) , nvecx, sc_d( 1, n_start ), nvecx, nbase, n_end-n_start+1 )

    array_out(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1) =  &
                array_in(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1)
#endif
    !
  END SUBROUTINE h2d_memsync_r2d
  !
  SUBROUTINE h2d_memsync_r3d(array_out, array_in, range1, range2, range3 )
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:,:)
    REAL(DP), INTENT(IN)    :: array_in(:,:,:)
    INTEGER, INTENT(IN) ::  range1(3), range2(3), range3(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: d2_start, d2_size, d2_ld
    INTEGER :: d3_start, d3_size, d3_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    d2_start = range2(1)
    d2_size = range2(2)
    d2_ld = range2(3)
    d3_start = range3(1)
    d3_size = range3(2)
    d3_ld = range3(3)
    !
#if defined(__CUDA)
    CALL errore('cu_memsync_','3D arrays not implemented yet',1)
#else
    CALL errore('cu_memsync_','3D arrays not implemented yet',1)
#endif
    !
  END SUBROUTINE h2d_memsync_r3d
  !
  SUBROUTINE h2d_memsync_c1d(array_out, array_in, range1 )
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:)
    COMPLEX(DP), INTENT(IN)    :: array_in(:)
    INTEGER, INTENT(IN) ::  range1(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyHostToDevice )
#else
    array_out(d1_start:d1_start+d1_size-1) = array_in(d1_start:d1_start+d1_size-1)
#endif
    !
  END SUBROUTINE h2d_memsync_c1d
  !
  SUBROUTINE h2d_memsync_c2d(array_out, array_in, range1, range2 )
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:)
    COMPLEX(DP), INTENT(IN)    :: array_in(:,:)
    INTEGER, INTENT(IN) ::  range1(3), range2(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: d2_start, d2_size, d2_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    d2_start = range2(1)
    d2_size = range2(2)
    d2_ld = range2(3)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else

    !pinned_buffer(1:nbase, n_start:n_end) = sc_d( 1:nbase, n_start:n_end )
    !ierr = cudaMemcpy2D( pinned_buffer(1, n_start) , nvecx, sc_d( 1, n_start ), nvecx, nbase, n_end-n_start+1 )

    array_out(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1) =  &
                array_in(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1)
#endif
    !
  END SUBROUTINE h2d_memsync_c2d
  !
  SUBROUTINE h2d_memsync_c3d(array_out, array_in, range1, range2, range3 )
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:,:)
    COMPLEX(DP), INTENT(IN)    :: array_in(:,:,:)
    INTEGER, INTENT(IN) ::  range1(3), range2(3), range3(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_out
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: d2_start, d2_size, d2_ld
    INTEGER :: d3_start, d3_size, d3_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    d2_start = range2(1)
    d2_size = range2(2)
    d2_ld = range2(3)
    d3_start = range3(1)
    d3_size = range3(2)
    d3_ld = range3(3)
    !
#if defined(__CUDA)
    CALL errore('cu_memsync_','3D arrays not implemented yet',1)
#else
    CALL errore('cu_memsync_','3D arrays not implemented yet',1)
#endif
    !
  END SUBROUTINE h2d_memsync_c3d
  !
!
  SUBROUTINE d2h_memsync_r1d(array_out, array_in, range1 )
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:)
    REAL(DP), INTENT(IN)    :: array_in(:)
    INTEGER, INTENT(IN) ::  range1(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_in
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyDeviceToHost )
#else
    array_out(d1_start:d1_start+d1_size-1) = array_in(d1_start:d1_start+d1_size-1)
#endif
    !
  END SUBROUTINE d2h_memsync_r1d
  !
  SUBROUTINE d2h_memsync_r2d(array_out, array_in, range1, range2 )
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:)
    REAL(DP), INTENT(IN)    :: array_in(:,:)
    INTEGER, INTENT(IN) ::  range1(3), range2(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_in
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: d2_start, d2_size, d2_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    d2_start = range2(1)
    d2_size = range2(2)
    d2_ld = range2(3)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else

    !pinned_buffer(1:nbase, n_start:n_end) = sc_d( 1:nbase, n_start:n_end )
    !ierr = cudaMemcpy2D( pinned_buffer(1, n_start) , nvecx, sc_d( 1, n_start ), nvecx, nbase, n_end-n_start+1 )

    array_out(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1) =  &
                array_in(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1)
#endif
    !
  END SUBROUTINE d2h_memsync_r2d
  !
  SUBROUTINE d2h_memsync_r3d(array_out, array_in, range1, range2, range3 )
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(INOUT) :: array_out(:,:,:)
    REAL(DP), INTENT(IN)    :: array_in(:,:,:)
    INTEGER, INTENT(IN) ::  range1(3), range2(3), range3(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_in
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: d2_start, d2_size, d2_ld
    INTEGER :: d3_start, d3_size, d3_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    d2_start = range2(1)
    d2_size = range2(2)
    d2_ld = range2(3)
    d3_start = range3(1)
    d3_size = range3(2)
    d3_ld = range3(3)
    !
#if defined(__CUDA)
    CALL errore('cu_memsync_','3D arrays not implemented yet',1)
#else
    CALL errore('cu_memsync_','3D arrays not implemented yet',1)
#endif
    !
  END SUBROUTINE d2h_memsync_r3d
  !
  SUBROUTINE d2h_memsync_c1d(array_out, array_in, range1 )
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:)
    COMPLEX(DP), INTENT(IN)    :: array_in(:)
    INTEGER, INTENT(IN) ::  range1(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_in
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyDeviceToHost )
#else
    array_out(d1_start:d1_start+d1_size-1) = array_in(d1_start:d1_start+d1_size-1)
#endif
    !
  END SUBROUTINE d2h_memsync_c1d
  !
  SUBROUTINE d2h_memsync_c2d(array_out, array_in, range1, range2 )
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:)
    COMPLEX(DP), INTENT(IN)    :: array_in(:,:)
    INTEGER, INTENT(IN) ::  range1(3), range2(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_in
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: d2_start, d2_size, d2_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    d2_start = range2(1)
    d2_size = range2(2)
    d2_ld = range2(3)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else

    !pinned_buffer(1:nbase, n_start:n_end) = sc_d( 1:nbase, n_start:n_end )
    !ierr = cudaMemcpy2D( pinned_buffer(1, n_start) , nvecx, sc_d( 1, n_start ), nvecx, nbase, n_end-n_start+1 )

    array_out(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1) =  &
                array_in(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1)
#endif
    !
  END SUBROUTINE d2h_memsync_c2d
  !
  SUBROUTINE d2h_memsync_c3d(array_out, array_in, range1, range2, range3 )
    !
#if defined(__CUDA)
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(INOUT) :: array_out(:,:,:)
    COMPLEX(DP), INTENT(IN)    :: array_in(:,:,:)
    INTEGER, INTENT(IN) ::  range1(3), range2(3), range3(3)
    !
#if defined(__CUDA)
    attributes(DEVICE) :: array_in
#endif
    !
    INTEGER :: d1_start, d1_size, d1_ld
    INTEGER :: d2_start, d2_size, d2_ld
    INTEGER :: d3_start, d3_size, d3_ld
    INTEGER :: ierr
    !
    d1_start = range1(1)
    d1_size = range1(2)
    d1_ld = range1(3)
    d2_start = range2(1)
    d2_size = range2(2)
    d2_ld = range2(3)
    d3_start = range3(1)
    d3_size = range3(2)
    d3_ld = range3(3)
    !
#if defined(__CUDA)
    CALL errore('cu_memsync_','3D arrays not implemented yet',1)
#else
    CALL errore('cu_memsync_','3D arrays not implemented yet',1)
#endif
    !
  END SUBROUTINE d2h_memsync_c3d
  !
END MODULE cuda_util
!
!{%- raw %}
!
!
! === TEMPLATE USED TO GENERATE THIS FILE ===
!
!    MODULE cuda_util
!      !
!      USE util_param,   ONLY : DP
!      !
!      IMPLICIT NONE
!      !
!      PUBLIC :: cuf_memcpy, cuf_memset, cu_memsync
!      !
!      INTERFACE cuf_memcpy
!        MODULE PROCEDURE &
!          {%- for t in types %}
!          {%- for d in range(dimensions) %}
!          cuf_memcpy_{{t[0]|lower}}{{d+1}}d{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}
!      END INTERFACE
!      !
!      INTERFACE cuf_memset
!        MODULE PROCEDURE &
!          {%- for t in types %}
!          {%- for d in range(dimensions) %}
!          cuf_memset_{{t[0]|lower}}{{d+1}}d{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}
!      END INTERFACE
!      !
!      INTERFACE cu_memsync
!        MODULE PROCEDURE &
!          {%- for t in types %}
!          {%- for d in range(dimensions) %}
!          h2d_memsync_{{t[0]|lower}}{{d+1}}d{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}
!    #if defined(__CUDA)
!        MODULE PROCEDURE &
!          {%- for t in types %}
!          {%- for d in range(dimensions) %}
!          d2h_memsync_{{t[0]|lower}}{{d+1}}d{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}{% if not loop.last %}, &{%- endif %}
!          {%- endfor %}
!    #endif
!      END INTERFACE
!      !
!      CONTAINS
!      !
!    {%- for t in types %}
!    {%- for d in range(1,dimensions+1) %}
!      SUBROUTINE cuf_memcpy_{{t[0]|lower}}{{d}}d(array_out, array_in,{% for dd in range(d) %} {{ "range%s"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %} )
!        !
!        IMPLICIT NONE
!        !
!        {{t}}(DP), INTENT(INOUT) :: array_out({% for dd in range(d) %}:{% if not loop.last %}, {%- endif %}{% endfor %})
!        {{t}}(DP), INTENT(IN) :: array_in({% for dd in range(d) %}:{% if not loop.last %}, {%- endif %}{% endfor %})
!        INTEGER, INTENT(IN) :: {% for dd in range(d) %} {{ "range%s(2)"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %}
!        !
!    #if defined(__CUDA)
!        attributes(DEVICE) :: array_out, array_in
!    #endif
!        !
!    {%- for dd in range(d) %}
!        INTEGER :: i{{dd+1}}, d{{dd+1}}s, d{{dd+1}}e
!    {%- endfor %}
!        !
!    {%- for dd in range(d) %}
!        d{{dd+1}}s = range{{dd+1}}(1)
!        d{{dd+1}}e = range{{dd+1}}(2)
!    {%- endfor %}
!        !
!        !$cuf kernel do({{d}})
!    {%- for dd in range(d,0,-1) %}
!        DO i{{dd}} = d{{dd}}s, d{{dd}}e
!    {%- endfor %}
!           array_out( {%- for dd in range(d) %}i{{dd+1}}{% if not loop.last %}, {%- endif %} {%- endfor %} ) = array_in( {%- for dd in range(d) %}i{{dd+1}}{% if not loop.last %}, {%- endif %} {%- endfor %} )
!    {%- for dd in range(d) %}
!        ENDDO
!    {%- endfor %}
!        !
!      END SUBROUTINE cuf_memcpy_{{t[0]|lower}}{{d}}d
!      !
!    {%- endfor %}
!    {%- endfor %}
!      !
!    {%- for t in types %}
!    {%- for d in range(1,dimensions+1) %}
!      SUBROUTINE cuf_memset_{{t[0]|lower}}{{d}}d(array_out, val,{% for dd in range(d) %} {{ "range%s"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %} )
!        !
!        IMPLICIT NONE
!        !
!        {{t}}(DP), INTENT(INOUT) :: array_out({% for dd in range(d) %}:{% if not loop.last %}, {%- endif %}{% endfor %})
!        {{t}}(DP), INTENT(IN) :: val
!        INTEGER, INTENT(IN) :: {% for dd in range(d) %} {{ "range%s(2)"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %}
!        !
!    #if defined(__CUDA)
!        attributes(DEVICE) :: array_out
!    #endif
!        !
!    {%- for dd in range(d) %}
!        INTEGER :: i{{dd+1}}, d{{dd+1}}s, d{{dd+1}}e
!    {%- endfor %}
!        !
!    {%- for dd in range(d) %}
!        d{{dd+1}}s = range{{dd+1}}(1)
!        d{{dd+1}}e = range{{dd+1}}(2)
!    {%- endfor %}
!        !
!        !$cuf kernel do({{d}})
!    {%- for dd in range(d,0,-1) %}
!        DO i{{dd}} = d{{dd}}s, d{{dd}}e
!    {%- endfor %}
!           array_out( {%- for dd in range(d) %}i{{dd+1}}{% if not loop.last %}, {%- endif %} {%- endfor %} ) = val
!    {%- for dd in range(d) %}
!        ENDDO
!    {%- endfor %}
!        !
!      END SUBROUTINE cuf_memset_{{t[0]|lower}}{{d}}d
!      !
!    {%- endfor %}
!    {%- endfor %}
!      !
!    {%- for t in types %}
!    {%- for d in range(1,dimensions+1) %}
!      SUBROUTINE h2d_memsync_{{t[0]|lower}}{{d}}d(array_out, array_in,{% for dd in range(d) %} {{ "range%s"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %} )
!    #if defined(__CUDA)
!        USE cudafor
!    #endif
!        !
!        IMPLICIT NONE
!        !
!        {{t}}(DP), INTENT(INOUT) :: array_out({% for dd in range(d) %}:{% if not loop.last %}, {%- endif %}{% endfor %})
!        {{t}}(DP), INTENT(IN)    :: array_in({% for dd in range(d) %}:{% if not loop.last %}, {%- endif %}{% endfor %})
!        INTEGER, INTENT(IN) :: {% for dd in range(d) %} {{ "range%s(3)"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %}
!        !
!    #if defined(__CUDA)
!        attributes(DEVICE) :: array_out
!    #endif
!        !
!    {%- for dd in range(d) %}
!        INTEGER :: d{{dd+1}}_start, d{{dd+1}}_size, d{{dd+1}}_ld
!    {%- endfor %}
!        INTEGER :: ierr
!        !
!    {%- for dd in range(d) %}
!        d{{dd+1}}_start = range{{dd+1}}(1)
!        d{{dd+1}}_size = range{{dd+1}}(2)
!        d{{dd+1}}_ld = range{{dd+1}}(3)
!    {%- endfor %}
!        !
!    #if defined(__CUDA)
!    {%- if d==1 %}
!        ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyHostToDevice )
!    {%- elif d==2 %}
!        ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
!    {%- elif d==3 %}
!        CALL errore('cu_memsync_','3D arrays not implemented yet',1)
!    {%- endif %}
!    #else
!    {%- if d==1 %}
!        array_out(d1_start:d1_start+d1_size-1) = array_in(d1_start:d1_start+d1_size-1)
!    {%- elif d==2 %}
!    
!        !pinned_buffer(1:nbase, n_start:n_end) = sc_d( 1:nbase, n_start:n_end )
!        !ierr = cudaMemcpy2D( pinned_buffer(1, n_start) , nvecx, sc_d( 1, n_start ), nvecx, nbase, n_end-n_start+1 )
!    
!        array_out(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1) =  &
!                    array_in(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1)
!    {%- elif d==3 %}
!        CALL errore('cu_memsync_','3D arrays not implemented yet',1)
!    {%- endif %}
!    #endif
!        !
!      END SUBROUTINE h2d_memsync_{{t[0]|lower}}{{d}}d
!      !
!    {%- endfor %}
!    {%- endfor %}
!    !
!    {%- for t in types %}
!    {%- for d in range(1,dimensions+1) %}
!      SUBROUTINE d2h_memsync_{{t[0]|lower}}{{d}}d(array_out, array_in,{% for dd in range(d) %} {{ "range%s"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %} )
!        !
!    #if defined(__CUDA)
!        USE cudafor
!    #endif
!        !
!        IMPLICIT NONE
!        !
!        {{t}}(DP), INTENT(INOUT) :: array_out({% for dd in range(d) %}:{% if not loop.last %}, {%- endif %}{% endfor %})
!        {{t}}(DP), INTENT(IN)    :: array_in({% for dd in range(d) %}:{% if not loop.last %}, {%- endif %}{% endfor %})
!        INTEGER, INTENT(IN) :: {% for dd in range(d) %} {{ "range%s(3)"|format(dd+1) }}{% if not loop.last %}, {%- endif %}{% endfor %}
!        !
!    #if defined(__CUDA)
!        attributes(DEVICE) :: array_in
!    #endif
!        !
!    {%- for dd in range(d) %}
!        INTEGER :: d{{dd+1}}_start, d{{dd+1}}_size, d{{dd+1}}_ld
!    {%- endfor %}
!        INTEGER :: ierr
!        !
!    {%- for dd in range(d) %}
!        d{{dd+1}}_start = range{{dd+1}}(1)
!        d{{dd+1}}_size = range{{dd+1}}(2)
!        d{{dd+1}}_ld = range{{dd+1}}(3)
!    {%- endfor %}
!        !
!    #if defined(__CUDA)
!    {%- if d==1 %}
!        ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyDeviceToHost )
!    {%- elif d==2 %}
!        ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
!    {%- elif d==3 %}
!        CALL errore('cu_memsync_','3D arrays not implemented yet',1)
!    {%- endif %}
!    #else
!    {%- if d==1 %}
!        array_out(d1_start:d1_start+d1_size-1) = array_in(d1_start:d1_start+d1_size-1)
!    {%- elif d==2 %}
!    
!        !pinned_buffer(1:nbase, n_start:n_end) = sc_d( 1:nbase, n_start:n_end )
!        !ierr = cudaMemcpy2D( pinned_buffer(1, n_start) , nvecx, sc_d( 1, n_start ), nvecx, nbase, n_end-n_start+1 )
!    
!        array_out(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1) =  &
!                    array_in(d1_start:d1_start+d1_size-1, d2_start:d2_start+d2_size-1)
!    {%- elif d==3 %}
!        CALL errore('cu_memsync_','3D arrays not implemented yet',1)
!    {%- endif %}
!    #endif
!        !
!      END SUBROUTINE d2h_memsync_{{t[0]|lower}}{{d}}d
!      !
!    {%- endfor %}
!    {%- endfor %}
!    END MODULE cuda_util
!
! === CODE TO GENERATE THE f90 FILE ===
!
!import sys, os, jinja2
!
!def render(tpl_path, context):
!    path, filename = os.path.split(tpl_path)
!    return jinja2.Environment(undefined=jinja2.StrictUndefined,
!        loader=jinja2.FileSystemLoader(path or './')
!    ).get_template(filename).render(context)
!with open('cuda_util.f90', 'w') as f: f.write(render('cuda_util.jf90', {'types': ['REAL', 'COMPLEX'], 'dimensions': 3}))
!
! {%- endraw %}
