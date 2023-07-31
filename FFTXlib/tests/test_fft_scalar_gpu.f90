#if defined(__CUDA)
#define devattr ,device
#else
#define devattr
#endif

#if defined(__CUDA) || defined(__OPENMP_GPU)
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
#if defined(__CUDA)
    USE cudafor
#endif
    USE fft_param, ONLY : DP
    implicit none
    complex(DP) devattr :: c_d(:)
    complex(DP)         :: c(:)
    integer, intent(in) :: n
    !
    real(DP), ALLOCATABLE :: rnd_aux(:)
    !
    ALLOCATE (rnd_aux(2*n))
    CALL RANDOM_NUMBER(rnd_aux)
    c = CMPLX(rnd_aux(1:n), rnd_aux(n:2*n))
    c_d = c
    !$omp target update to(c_d)
    DEALLOCATE(rnd_aux)
  END SUBROUTINE fill_random
  !
  SUBROUTINE test_cft_1z_gpu(test)
    USE fft_scalar, ONLY : cft_1z, cft_2xy, cfft3d, cfft3ds
    USE fft_param, ONLY : DP
#if defined(__CUDA)
    USE cudafor
    USE fft_scalar, ONLY : cft_1z_gpu, cft_2xy_gpu, cfft3d_gpu, cfft3ds_gpu
#elif defined(__OPENMP_GPU)
    USE fft_scalar, ONLY : cft_1z_omp, cft_2xy_omp, cfft3d_omp, cfft3ds_omp
#endif
    implicit none
    !
    TYPE(tester_t) :: test
    !
#if !defined(__OPENMP_GPU)
    ! stream
    integer(kind = cuda_stream_kind) :: stream = 0
#endif
    ! size
    integer, parameter :: nsl=100, nz=56, ldz=56
    !
    ! array for random values
    real(DP) :: rnd_aux(2 * nsl * ldz)
    !
    ! variables on device
    complex(DP) devattr :: c_d(nsl*ldz)
    complex(DP) devattr :: cout_d(nsl*ldz)
    ! variables on host
    complex(DP) :: c(nsl*ldz)
    complex(DP) :: cout(nsl*ldz)
    !
    !$omp target enter data map(alloc:c_d)
    !$omp target enter data map(alloc:cout_d)
    !
    CALL fill_random(c, c_d, nsl*ldz)
    !
    ! Check forward direction
    CALL cft_1z(c, nsl, nz, ldz, 1, cout)
#if !defined(__OPENMP_GPU)
    CALL cft_1z_gpu(c_d, nsl, nz, ldz, 1, cout_d, stream)
#else
    CALL cft_1z_omp(c_d, nsl, nz, ldz, 1, cout_d)
#endif
    !$omp target update from(cout_d)
    !
#if !defined(__OPENMP_GPU)
    ! Use c as auxiliary variable hosting GPU results
    c = (0.d0, 0.d0)
    c = cout_d
    CALL test%assert_close( c, cout )
#else
    CALL test%assert_close( cout_d, cout )
#endif
    !
    !
    CALL fill_random(c, c_d, nsl*ldz)
    !
    ! Check backward direction
    CALL cft_1z(c, nsl, nz, ldz, -1, cout)
#if !defined(__OPENMP_GPU)
    CALL cft_1z_gpu(c_d, nsl, nz, ldz, -1, cout_d, stream)
#else
    CALL cft_1z_omp(c_d, nsl, nz, ldz, -1, cout_d)
#endif
    !$omp target exit data map(from:cout_d)
    !
#if !defined(__OPENMP_GPU)
    ! Use c as auxiliary variable hosting GPU results
    c = (0.d0, 0.d0)
    c = cout_d
    CALL test%assert_close( c, cout )
#else
    CALL test%assert_close( cout_d, cout )
#endif
    !
    ! ==  Same as above, for inplace call ==
    !
    CALL fill_random(c, c_d, nsl*ldz)
    !
    ! Check forward direction
    CALL cft_1z(c, nsl, nz, ldz, 1, cout)
#if !defined(__OPENMP_GPU)
    CALL cft_1z_gpu(c_d, nsl, nz, ldz, 1, cout_d, stream, in_place=.true.)
#else
    CALL cft_1z_omp(c_d, nsl, nz, ldz, 1, cout_d, in_place=.true.)
#endif
    !$omp target update from(c_d)
    !
#if !defined(__OPENMP_GPU)
    ! Use c as auxiliary variable hosting GPU results
    c = (0.d0, 0.d0)
    c = c_d
    CALL test%assert_close( c, cout )
#else
    CALL test%assert_close( c_d, cout )
#endif
    !
    !
    CALL fill_random(c, c_d, nsl*ldz)
    !
    ! Check backward direction
    CALL cft_1z(c, nsl, nz, ldz, -1, cout)
#if !defined(__OPENMP_GPU)
    CALL cft_1z_gpu(c_d, nsl, nz, ldz, -1, cout_d, stream, in_place=.true.)
#else
    CALL cft_1z_omp(c_d, nsl, nz, ldz, -1, cout_d, in_place=.true.)
#endif
    !$omp target exit data map(from:c_d)
    !
#if !defined(__OPENMP_GPU)
    ! Use c as auxiliary variable hosting GPU results
    c = (0.d0, 0.d0)
    c = c_d
    CALL test%assert_close( c, cout )
#else
    CALL test%assert_close( c_d, cout )
#endif
    !
  END SUBROUTINE test_cft_1z_gpu
  !
  SUBROUTINE test_cft_2xy_gpu(test)
    USE fft_scalar, ONLY : cft_2xy
    USE fft_param, ONLY : DP
#if defined(__CUDA)
    USE cudafor
    USE fft_scalar, ONLY : cft_2xy_gpu
#elif defined(__OPENMP_GPU)
    USE fft_scalar, ONLY : cft_2xy_omp
#endif
    implicit none
    !
    TYPE(tester_t) :: test
    !
#if !defined(__OPENMP_GPU)
    ! stream
    integer(kind = cuda_stream_kind) :: stream = 0
#endif
    ! size
    integer, parameter :: nx = 10, ny = 10, nzl = 5, ldx=10, ldy=10
    !
    ! array for random values
    real(DP) :: rnd_aux(2 * nzl * ldx * ldy)
    ! variables on device
    complex(DP) devattr :: c_d(nzl * ldx * ldy)
    complex(DP) devattr :: tmp_d(nzl * ldx * ldy)
    ! variables on host
    complex(DP) :: c(nzl * ldx * ldy)
    complex(DP) :: tmp(nzl * ldx * ldy)
    !
    !$omp target enter data map(alloc:c_d)
    !
    CALL fill_random(c, c_d, nzl * ldx * ldy)
    !
#if !defined(__OPENMP_GPU)
    CALL cft_2xy_gpu(c_d, tmp_d, nzl, nx, ny, ldx, ldy, 1, stream)
#else
    CALL cft_2xy_omp(c_d, nzl, nx, ny, ldx, ldy, 1)
#endif
    !$omp target update from(c_d)
    CALL cft_2xy(c, nzl, nx, ny, ldx, ldy, 1)
    !
#if !defined(__OPENMP_GPU)
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
    CALL test%assert_close( c, tmp )
#else
    CALL test%assert_close( c, c_d )
#endif
    !
    CALL fill_random(c, c_d, nzl * ldx * ldy)
    !
#if !defined(__OPENMP_GPU)
    CALL cft_2xy_gpu(c_d, tmp_d, nzl, nx, ny, ldx, ldy, -1, stream)
#else
    CALL cft_2xy_omp(c_d, nzl, nx, ny, ldx, ldy, -1)
#endif
    !$omp target exit data map(from:c_d)
    CALL cft_2xy(c, nzl, nx, ny, ldx, ldy, -1)
    !
#if !defined(__OPENMP_GPU)
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
    CALL test%assert_close( c, tmp )
#else
    CALL test%assert_close( c, c_d )

#endif
    !
  END SUBROUTINE test_cft_2xy_gpu
  !
  SUBROUTINE test_cfft3d_gpu(test)
    USE fft_scalar, ONLY : cfft3d
    USE fft_param, ONLY : DP
#if defined(__CUDA)
    USE cudafor
    USE fft_scalar, ONLY : cfft3d_gpu
#elif defined(__OPENMP_GPU)
    USE fft_scalar, ONLY : cfft3d_omp
#endif
    implicit none
    !
    TYPE(tester_t) :: test
    !
#if !defined(__OPENMP_GPU)
    ! stream
    integer(kind = cuda_stream_kind) :: stream = 0
#endif
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
    complex(DP) devattr :: c_d(howmany * ldx * ldy * ldz)
    complex(DP) devattr :: tmp_d(howmany * ldx * ldy * ldz)
    ! variables on host
    complex(DP) :: c(howmany * ldx * ldy * ldz)
    complex(DP) :: tmp(howmany * ldx * ldy * ldz)
    !
    integer :: i, rs, re
    !
    !$omp target enter data map(alloc:c_d)
    !
#if ! defined(__DFTI)
    print *, 'The current CPU scalar driver does not support howmany. Reverting to howmany 1'
#endif

    CALL fill_random(c, c_d, howmany * ldx * ldy * ldz)
    !
    CALL cfft3d( c, nx, ny, nz, ldx, ldy, ldz, howmany, 1 )
#if !defined(__OPENMP_GPU)
    CALL cfft3d_gpu( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, 1, stream )
#else
    CALL cfft3d_omp( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, 1 )
#endif
    !$omp target update from(c_d)
    !
#if !defined(__OPENMP_GPU)
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
#endif
    DO i=0, howmany - 1
      rs = i * ldx * ldy * ldz + 1
      re = rs + nx * ny  * nz  - 1
#if !defined(__OPENMP_GPU)
      CALL test%assert_close( c(rs:re), tmp(rs:re) )
#else
      CALL test%assert_close( c(rs:re), c_d(rs:re) )
#endif
    END DO
    !
    CALL fill_random(c, c_d, howmany * ldx * ldy * ldz)
    !
    CALL cfft3d( c, nx, ny, nz, ldx, ldy, ldz, howmany, -1 )
#if !defined(__OPENMP_GPU)
    CALL cfft3d_gpu( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, -1, stream )
#else
    CALL cfft3d_omp( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, -1 )
#endif
    !$omp target exit data map(from:c_d)
    !
#if !defined(__OPENMP_GPU)
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
#endif
    DO i=0, howmany - 1
      rs = i * ldx * ldy * ldz + 1
      re = rs + nx * ny  * nz  - 1
#if !defined(__OPENMP_GPU)
      CALL test%assert_close( c(rs:re), tmp(rs:re) )
#else
      CALL test%assert_close( c(rs:re), c_d(rs:re) )
#endif
    END DO
    !
  END SUBROUTINE test_cfft3d_gpu
  !
  SUBROUTINE test_cfft3ds_gpu(test)
    USE fft_scalar, ONLY : cfft3ds
    USE fft_param, ONLY : DP
#if defined(__CUDA)
    USE cudafor
    USE fft_scalar, ONLY : cfft3ds_gpu
#elif defined(__OPENMP_GPU)
    USE fft_scalar, ONLY : cfft3ds_omp
#endif
    implicit none
    !
    TYPE(tester_t) :: test
    !
#if !defined(__OPENMP_GPU)
    ! stream
    integer(kind = cuda_stream_kind) :: stream = 0
#endif
    ! size
    integer, parameter :: nx = 10, ny = 10, nz = 10
    integer, parameter :: ldx= 10, ldy= 10, ldz= 10, howmany=1
    !
    ! array for random values
    real(DP) :: rnd_aux(2 * howmany * ldx * ldy * ldz)
    ! variables on device
    complex(DP) devattr :: c_d(howmany * ldx * ldy * ldz)
    complex(DP) devattr :: tmp_d(howmany * ldx * ldy * ldz)
    ! variables on host
    complex(DP) :: c(howmany * ldx * ldy * ldz)
    complex(DP) :: tmp(howmany * ldx * ldy * ldz)
    integer     :: do_fft_y(ldx), do_fft_z(ldx*ldy)
    !
    !$omp target enter data map(alloc:c_d)
    !
    CALL fill_random(c, c_d, howmany * ldx * ldy * ldz)
    do_fft_y = 1; do_fft_z = 1
    !
    CALL cfft3ds( c, nx, ny, nz, ldx, ldy, ldz, howmany, 1, do_fft_z, do_fft_y)
#if !defined(__OPENMP_GPU)
    CALL cfft3ds_gpu( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, 1, do_fft_z, do_fft_y, stream)
#else
    CALL cfft3ds_omp( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, 1, do_fft_z, do_fft_y)
#endif
    !$omp target update from(c_d)
    !
#if !defined(__OPENMP_GPU)
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
    CALL test%assert_close( c, tmp )
#else
    CALL test%assert_close( c, c_d )
#endif
    !
    CALL fill_random(c, c_d, howmany * ldx * ldy * ldz)
    !
    CALL cfft3ds( c, nx, ny, nz, ldx, ldy, ldz, howmany, -1, do_fft_z, do_fft_y)
#if !defined(__OPENMP_GPU)
    CALL cfft3ds_gpu( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, -1, do_fft_z, do_fft_y, stream )
#else
    CALL cfft3ds_omp( c_d, nx, ny, nz, ldx, ldy, ldz, howmany, -1, do_fft_z, do_fft_y )
#endif
    !$omp target exit data map(from:c_d)
    !
#if !defined(__OPENMP_GPU)
    ! Use c as auxiliary variable hosting GPU results
    tmp = c_d
    CALL test%assert_close( c, tmp )
#else
    CALL test%assert_close( c, c_d )
#endif
    !
  END SUBROUTINE test_cfft3ds_gpu

end program test_fft_scalar_gpu
#else
program test_fft_scalar_gpu
end program test_fft_scalar_gpu
#endif
!
! Dummy
SUBROUTINE stop_clock(label)
CHARACTER(*) :: label
END SUBROUTINE stop_clock
!
SUBROUTINE start_clock(label)
CHARACTER(*) :: label
END SUBROUTINE start_clock
!
SUBROUTINE stop_clock_gpu(label)
CHARACTER(*) :: label
END SUBROUTINE stop_clock_gpu
!
SUBROUTINE start_clock_gpu(label)
CHARACTER(*) :: label
END SUBROUTINE start_clock_gpu
