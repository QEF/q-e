#if defined(__CUDA)
program test_fft_scatter_mod_gpu
#if defined(__MPI)
    USE mpi
#endif
    USE tester
    IMPLICIT NONE
    ! MPI type
    type mpi_t
    integer :: me, n, root, comm
    end type mpi_t
    TYPE(mpi_t) :: mp
    !
    TYPE(tester_t) :: test
    !
    INTEGER :: ierr, level, i
    !
#if defined(__MPI)
#if defined(_OPENMP)
  CALL MPI_Init_thread(MPI_THREAD_FUNNELED,level, ierr)
#else
  CALL MPI_Init(ierr)
#endif
#endif
    !
    CALL mpi_data_init(mp%me, mp%n, mp%root, mp%comm)
    !
    CALL test%init()
    !
    test%tolerance64 = 1.d-14
    !
    CALL save_random_seed("test_fft_scatter_mod_gpu", mp%me)
    !
    DO i = 1, mp%n
      IF (MOD(mp%n,i) == 0 ) THEN
        CALL test_fft_scatter_xy_gpu_1(mp, test, .true., i)
        CALL test_fft_scatter_xy_gpu_1(mp, test, .false., i)
        !
        CALL test_fft_scatter_yz_gpu_1(mp, test, .true., i)
        CALL test_fft_scatter_yz_gpu_1(mp, test, .false., i)
      END IF
    END DO
    !
    CALL collect_results(test)
    !
    IF (mp%me == mp%root) CALL test%print()
    !
#if defined(__MPI)
    CALL MPI_Finalize(ierr)
#endif
  CONTAINS
  !
  SUBROUTINE mpi_data_init(mpme, npes, mproot, comm)
    implicit none
    integer, intent(out) :: mpme, npes, mproot, comm
    integer :: ierr
    mpme=0; npes=1; mproot=0; comm=0
#if defined(__MPI)
    CALL mpi_comm_rank(MPI_COMM_WORLD, mpme, ierr)
    CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)
    comm = MPI_COMM_WORLD
#endif
  END SUBROUTINE mpi_data_init
  SUBROUTINE fft_desc_init(dfft, smap, flavor, gamma_only, parallel, comm, nyfft)
    USE fft_types,       ONLY : fft_type_descriptor, fft_type_init
    USE stick_base,      ONLY : sticks_map
    USE fft_param, ONLY : DP
    implicit none
    TYPE(fft_type_descriptor) :: dfft
    TYPE(sticks_map) :: smap
    CHARACTER(LEN=*), INTENT(IN) :: flavor
    LOGICAL :: gamma_only
    LOGICAL :: parallel
    INTEGER :: comm, nyfft
    REAL(DP), PARAMETER :: pi=4.D0*DATAN(1.D0)
    !
    REAL(DP) :: at(3:3), bg(3:3)
    !
    at = RESHAPE((/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/), shape(at))
    bg = RESHAPE((/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/), shape(bg))
    bg = 2.d0*pi
    !
    CALL fft_type_init(dfft, smap, flavor, gamma_only, parallel, comm, at, bg, 12.d0, 6.d0, nyfft=nyfft)
    !
  END SUBROUTINE fft_desc_init
  
  SUBROUTINE fft_desc_finalize(dfft, smap)
    USE fft_types,       ONLY : fft_type_descriptor, fft_type_deallocate
    USE stick_base,      ONLY : sticks_map, sticks_map_deallocate
    implicit none
    TYPE(fft_type_descriptor) :: dfft
    TYPE(sticks_map) :: smap
    !
    CALL fft_type_deallocate(dfft)
    CALL sticks_map_deallocate( smap )
  END SUBROUTINE fft_desc_finalize
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
  SUBROUTINE test_fft_scatter_xy_gpu_1(mp, test, gamma_only, ny)
    USE cudafor
    USE fft_param,       ONLY : DP
    USE fft_types,       ONLY : fft_type_descriptor
    USE stick_base,      ONLY : sticks_map
    USE scatter_mod,     ONLY : fft_scatter_xy
    USE scatter_mod_gpu, ONLY : fft_scatter_xy_gpu
    implicit none
    TYPE(mpi_t) :: mp
    TYPE(tester_t) :: test
    !
    TYPE(fft_type_descriptor) :: dfft
    TYPE(sticks_map) :: smap
    LOGICAL, INTENT(IN) :: gamma_only
    INTEGER, INTENT(IN) :: ny
    !
    LOGICAL :: parallel
    COMPLEX(DP), ALLOCATABLE :: scatter_in(:), scatter_out(:), aux(:)
    COMPLEX(DP), ALLOCATABLE, DEVICE :: scatter_in_d(:), scatter_out_d(:)
    integer(kind = cuda_stream_kind) :: stream = 0
    !
    parallel = mp%n .gt. 1
    CALL fft_desc_init(dfft, smap, "wave", gamma_only, parallel, mp%comm, nyfft=ny)
    !
    ! Allocate variables
    ALLOCATE(scatter_in(dfft%nnr), scatter_out(dfft%nnr), aux(dfft%nnr))
    ALLOCATE(scatter_in_d(dfft%nnr), scatter_out_d(dfft%nnr))
    !
    ! Test 1
    CALL fill_random(scatter_in, scatter_in_d, dfft%nnr)
    !
    CALL fft_scatter_xy( dfft, scatter_in, scatter_out, dfft%nnr, 1 )
    CALL fft_scatter_xy_gpu( dfft, scatter_in_d, scatter_out_d, dfft%nnr, 1, stream )
    aux = scatter_out_d
    ! Check
    CALL test%assert_close( scatter_out, aux )
    !
    ! Test 2
    CALL fill_random(scatter_in, scatter_in_d, dfft%nnr)
    !
    CALL fft_scatter_xy( dfft, scatter_out, scatter_in, dfft%nnr, -1 )
    CALL fft_scatter_xy_gpu( dfft, scatter_out_d, scatter_in_d, dfft%nnr, -1, stream )
    aux = scatter_out_d
    ! Check
    CALL test%assert_close( scatter_out, aux )
    !
    CALL fft_desc_finalize(dfft, smap)
    DEALLOCATE(scatter_in, scatter_out, aux, scatter_in_d, scatter_out_d)
    !
  END SUBROUTINE test_fft_scatter_xy_gpu_1
  !
  SUBROUTINE test_fft_scatter_yz_gpu_1(mp, test, gamma_only, ny)
    !
    ! This test checks wave fft scatter, with parallel = .true. if 
    !  called with more than 1 MPI.
    !
    USE cudafor
    USE fft_param,       ONLY : DP
    USE fft_types,       ONLY : fft_type_descriptor
    USE stick_base,      ONLY : sticks_map
    USE scatter_mod,     ONLY : fft_scatter_yz
    USE scatter_mod_gpu, ONLY : fft_scatter_yz_gpu
    implicit none
    TYPE(mpi_t) :: mp
    TYPE(tester_t) :: test
    !
    TYPE(fft_type_descriptor) :: dfft
    TYPE(sticks_map) :: smap
    LOGICAL, INTENT(IN) :: gamma_only
    INTEGER, INTENT(IN) :: ny
    !
    LOGICAL :: parallel
    COMPLEX(DP), ALLOCATABLE :: scatter_in(:), scatter_out(:), aux(:)
    COMPLEX(DP), ALLOCATABLE, DEVICE :: scatter_in_d(:), scatter_out_d(:)
    integer(kind = cuda_stream_kind) :: stream = 0
    !
    parallel = mp%n .gt. 1
    CALL fft_desc_init(dfft, smap, "wave", gamma_only, parallel, mp%comm, nyfft=ny)
    !
    ! Allocate variables
    ALLOCATE(scatter_in(dfft%nnr), scatter_out(dfft%nnr), aux(dfft%nnr))
    ALLOCATE(scatter_in_d(dfft%nnr), scatter_out_d(dfft%nnr))
    !
    ! Test 1
    CALL fill_random(scatter_in, scatter_in_d, dfft%nnr)
    !
    CALL fft_scatter_yz( dfft, scatter_in, scatter_out, dfft%nnr, 1 )
    CALL fft_scatter_yz_gpu( dfft, scatter_in_d, scatter_out_d, dfft%nnr, 1, stream )
    aux = scatter_out_d
    ! Check
    CALL test%assert_close( scatter_out, aux )
    !
    ! Test 2
    CALL fill_random(scatter_in, scatter_in_d, dfft%nnr)
    !
    CALL fft_scatter_yz( dfft, scatter_out, scatter_in, dfft%nnr, -1 )
    CALL fft_scatter_yz_gpu( dfft, scatter_out_d, scatter_in_d, dfft%nnr, -1, stream )
    aux = scatter_out_d
    ! Check
    CALL test%assert_close( scatter_out, aux )
    !
    CALL fft_desc_finalize(dfft, smap)
    DEALLOCATE(scatter_in, scatter_out, aux, scatter_in_d, scatter_out_d)
    !
  END SUBROUTINE test_fft_scatter_yz_gpu_1
  
end program test_fft_scatter_mod_gpu
! dummy subroutines
subroutine stop_clock( label )
character(len=*) :: label
end subroutine stop_clock
subroutine start_clock( label )
character(len=*) :: label
end subroutine start_clock
!
#else
program test_fft_scatter_mod_gpu
end program test_fft_scatter_mod_gpu
#endif
