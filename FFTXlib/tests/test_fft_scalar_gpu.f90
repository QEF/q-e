#if defined(__CUDA)
program test_fft_scalar_gpu
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    !
    CALL test%init()
    !
    test%tolerance64 = 1.d-13
    !
    CALL save_random_seed("test_fft_scalar_gpu", 0)
    !
    CALL test_cft_2xy_gpu(test)
    !
    CALL test_cft_1z_gpu(test)
    !
    CALL test_cfft3d_gpu(test)
    !
    CALL test_cfft3ds_gpu(test)
    !
    CALL test%print()
    !
  CONTAINS
  !
  SUBROUTINE fill_random(c, c_d, n)
    USE cudafor
    USE fft_param, ONLY : DP
    implicit none
    complex(DP), device :: c_d(:)
    complex(DP)         :: c(:)
    integer, intent(in) :: n
    !
    real(DP), ALLOCATABLE :: rnd_aux(:)
    !
    ALLOCATE (rnd_aux(2*n))
    CALL RANDOM_NUMBER(rnd_aux)
    c = CMPLX(rnd_aux(1:n), rnd_aux(n:2*n))
    c_d = c
    DEALLOCATE(rnd_aux)
  END SUBROUTINE fill_random
  !
  SUBROUTINE test_cft_1z_gpu(test)
    USE fft_scalar, ONLY : cft_1z, cft_2xy, cfft3d, cfft3ds
    USE fft_scalar, ONLY : cft_1z_gpu, cft_2xy_gpu, cfft3d_gpu, cfft3ds_gpu
    USE fft_param, ONLY : DP
    USE cudafor
    implicit none
    !
    TYPE(tester_t) :: test
    !
    ! stream
    integer(kind = cuda_stream_kind) :: stream = 0
    ! size
    integer, parameter :: nsl=100, nz=56, ldz=56
    !
    ! array for random values
    real(DP) :: rnd_aux(2 * nsl * ldz)
    !
    ! variables on device
    complex(DP), device :: c_d(nsl*ldz)
    complex(DP), device :: cout_d(nsl*ldz)
    ! variables on host
    complex(DP) :: c(nsl*ldz)
    complex(DP) :: cout(nsl*ldz)
    !
    CALL fill_random(c, c_d, nsl*ldz)
    !
    ! Check forward direction
    CALL cft_1z(c, nsl, nz, ldz, 1, cout)
    CALL cft_1z_gpu(c_d, nsl, nz, ldz, 1, cout_d, stream)
    !
    ! Use c as auxiliary variable hosting GPU results
    c = (0.d0, 0.d0)
    c = cout_d
    CALL test%assert_close( c, cout )
    !
    !
    CALL fill_random(c, c_d, nsl*ldz)
    !
    ! Check backward direction
    CALL cft_1z(c, nsl, nz, ldz, -1, cout)
    CALL cft_1z_gpu(c_d, nsl, nz, ldz, -1, cout_d, stream)
    !
    ! Use c as auxiliary variable hosting GPU results
    c = (0.d0, 0.d0)
    c = cout_d
    CALL test%assert_close( c, cout )
    !
    ! ==  Same as above, for inplace call ==
    !
    CALL fill_random(c, c_d, nsl*ldz)
    !
    ! Check forward direction
    CALL cft_1z(c, nsl, nz, ldz, 1, cout)
    CALL cft_1z_gpu(c_d, nsl, nz, ldz, 1, cout_d, stream, in_place=.true.)
    !
    ! Use c as auxiliary variable hosting GPU results
    c = (0.d0, 0.d0)
    c = c_d
    CALL test%assert_close( c, cout )
    !
    !
    CALL fill_random(c, c_d, nsl*ldz)
    !
    ! Check backward direction
    CALL cft_1z(c, nsl, nz, ldz, -1, cout)
    CALL cft_1z_gpu(c_d, nsl, nz, ldz, -1, cout_d, stream, in_place=.true.)
    !
    ! Use c as auxiliary variable hosting GPU results
    c = (0.d0, 0.d0)
    c = c_d
    CALL test%assert_close( c, cout )
    !
  END SUBROUTINE test_cft_1z_gpu
  !
  SUBROUTINE test_cft_2xy_gpu(test)
    USE fft_scalar, ONLY : cft_2xy
    USE fft_scalar, ONLY : cft_2xy_gpu
    USE fft_param, ONLY : DP
    USE cudafor
    implicit none
    !
    TYPE(tester_t) :: test
    !
    ! stream
    integer(kind = cuda_stream_kind) :: stream = 0
    ! size
    integer, parameter :: nx = 10, ny = 10, nzl = 5, ldx=10, ldy=10
    !
    ! array for random values
    real(DP) :: rnd_aux(2 * nzl * ldx * ldy)
    ! variables on device
    complex(DP), device :: c_d(nzl * ldx * ldy)
    complex(DP), device :: tmp_d(nzl * ldx * ldy)
    ! variables on host
    complex(DP) :: c(nzl * ldx * ldy)
    complex(DP) :: tmp(nzl * ldx * ldy)
    !
    CALL fill_random(c, c_d, nzl * ldx * ldy)
    !
    CALL cft_2xy_gpu(c_d, tmp_d, nzl, nx, ny, ldx, ldy, 1, stream)
    CALL cft_2xy(c, nzl, nx, ny, ldx, ldy, 1)
    !
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
    CALL test%assert_close( c, tmp )
    !
    CALL fill_random(c, c_d, nzl * ldx * ldy)
    !
    CALL cft_2xy_gpu(c_d, tmp_d, nzl, nx, ny, ldx, ldy, -1, stream)
    CALL cft_2xy(c, nzl, nx, ny, ldx, ldy, -1)
    !
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
    CALL test%assert_close( c, tmp )
    !
  END SUBROUTINE test_cft_2xy_gpu
  !
  SUBROUTINE test_cfft3d_gpu(test)
    USE fft_scalar, ONLY : cfft3d
    USE fft_scalar, ONLY : cfft3d_gpu
    USE fft_param, ONLY : DP
    USE cudafor
    implicit none
    !
    TYPE(tester_t) :: test
    !
    ! stream
    integer(kind = cuda_stream_kind) :: stream = 0
    ! size
    integer, parameter :: nx = 10, ny = 10, nz = 10
    integer, parameter :: ldx= 10, ldy= 10, ldz= 10
#if ! defined(__DFTI)
    integer, parameter :: howmany=1
#else
    integer, parameter :: howmany=12
#endif
    !
    ! array for random values
    real(DP) :: rnd_aux(2 * howmany * ldx * ldy * ldz)
    ! variables on device
    complex(DP), device :: c_d(howmany * ldx * ldy * ldz)
    complex(DP), device :: tmp_d(howmany * ldx * ldy * ldz)
    ! variables on host
    complex(DP) :: c(howmany * ldx * ldy * ldz)
    complex(DP) :: tmp(howmany * ldx * ldy * ldz)
    !
#if ! defined(__DFTI)
    print *, 'The current CPU scalar driver does not support howmany. Reverting to howmany 1'
#endif
    
    CALL fill_random(c, c_d, howmany * ldx * ldy * ldz)
    !
    CALL cfft3d( c, nx, ny, nz, ldx, ldy, ldz, howmany, 1 )
    CALL cfft3d_gpu( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, 1, stream )
    !
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
    CALL test%assert_close( c, tmp )
    !
    CALL fill_random(c, c_d, howmany * ldx * ldy * ldz)
    !
    CALL cfft3d( c, nx, ny, nz, ldx, ldy, ldz, howmany, -1 )
    CALL cfft3d_gpu( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, -1, stream )
    !
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
    CALL test%assert_close( c, tmp )
    !
  END SUBROUTINE test_cfft3d_gpu
  !
  SUBROUTINE test_cfft3ds_gpu(test)
    USE fft_scalar, ONLY : cfft3ds
    USE fft_scalar, ONLY : cfft3ds_gpu
    USE fft_param, ONLY : DP
    USE cudafor
    implicit none
    !
    TYPE(tester_t) :: test
    !
    ! stream
    integer(kind = cuda_stream_kind) :: stream = 0
    ! size
    integer, parameter :: nx = 10, ny = 10, nz = 10
    integer, parameter :: ldx= 10, ldy= 10, ldz= 10, howmany=1
    !
    ! array for random values
    real(DP) :: rnd_aux(2 * howmany * ldx * ldy * ldz)
    ! variables on device
    complex(DP), device :: c_d(howmany * ldx * ldy * ldz)
    complex(DP), device :: tmp_d(howmany * ldx * ldy * ldz)
    ! variables on host
    complex(DP) :: c(howmany * ldx * ldy * ldz)
    complex(DP) :: tmp(howmany * ldx * ldy * ldz)
    integer     :: do_fft_y(ldx), do_fft_z(ldx*ldy)
    !
    CALL fill_random(c, c_d, howmany * ldx * ldy * ldz)
    do_fft_y = 1; do_fft_z = 1
    !
    CALL cfft3ds( c, nx, ny, nz, ldx, ldy, ldz, howmany, 1, do_fft_z, do_fft_y)
    CALL cfft3ds_gpu( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, 1, do_fft_z, do_fft_y, stream)
    !
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
    CALL test%assert_close( c, tmp )
    !
    CALL fill_random(c, c_d, howmany * ldx * ldy * ldz)
    !
    CALL cfft3ds( c, nx, ny, nz, ldx, ldy, ldz, howmany, -1, do_fft_z, do_fft_y)
    CALL cfft3ds_gpu( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, -1, do_fft_z, do_fft_y, stream )
    !
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
    CALL test%assert_close( c, tmp )
    !
  END SUBROUTINE test_cfft3ds_gpu
  
end program test_fft_scalar_gpu
#else
program test_fft_scalar_gpu
end program test_fft_scalar_gpu
#endif
