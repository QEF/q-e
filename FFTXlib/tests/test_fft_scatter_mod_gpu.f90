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
        ! gamma case
        CALL test_fft_scatter_xy_gpu_1(mp, test, .true., i)
        ! k case
        CALL test_fft_scatter_xy_gpu_1(mp, test, .false., i)
        !
        ! gamma case
        CALL test_fft_scatter_yz_gpu_1(mp, test, .true., i)
        ! k case
        CALL test_fft_scatter_yz_gpu_1(mp, test, .false., i)
      END IF
    END DO
    CALL test_fft_scatter_many_yz_gpu_1(mp, test, .true., 1)
    CALL test_fft_scatter_many_yz_gpu_1(mp, test, .false., 1)

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
    REAL(DP) :: at(3,3), bg(3,3)
    !
    at = RESHAPE((/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/), shape(at))
    bg = RESHAPE((/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/), shape(bg))
    bg = 2.d0*pi
    !
    CALL fft_type_init(dfft, smap, flavor, gamma_only, parallel, comm, at, bg, 12.d0, 4.d0, nyfft=nyfft)
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
    integer :: fft_sign = 2
    integer :: vsiz, nr1p_, compare_len, me2
    !
    parallel = mp%n .gt. 1
    
    CALL fft_desc_init(dfft, smap, "wave", gamma_only, parallel, mp%comm, nyfft=ny)
    me2    = dfft%mype2 + 1
    vsiz = dfft%nnr
    compare_len = dfft%nr1x * dfft%my_nr2p * dfft%my_nr3p
    if (ny > 1) then
       ! When using task groups, wave FFTs are not distributed along Y
       fft_sign = 3
       vsiz = dfft%nnr_tg
       compare_len = dfft%nr1x * dfft%nr2x * dfft%my_nr3p
    end if
    !
    ! Allocate variables
    ALLOCATE(scatter_in(vsiz), scatter_out(vsiz), aux(vsiz))
    ALLOCATE(scatter_in_d(vsiz), scatter_out_d(vsiz))
    !
    ! Test 1
    CALL fill_random(scatter_in, scatter_in_d, vsiz)
    !
    CALL fft_scatter_xy( dfft, scatter_in, scatter_out, vsiz, fft_sign )
    CALL fft_scatter_xy_gpu( dfft, scatter_in_d, scatter_out_d, vsiz, fft_sign, stream )
    aux(1:compare_len) = scatter_out_d(1:compare_len)
    !
    ! Check
    CALL test%assert_close( scatter_out(1:compare_len), aux(1:compare_len) )
    !
    ! Test 2
    CALL fill_random(scatter_in, scatter_in_d, vsiz)
    !
    CALL fft_scatter_xy( dfft, scatter_out, scatter_in, vsiz, -1*fft_sign )
    CALL fft_scatter_xy_gpu( dfft, scatter_out_d, scatter_in_d, vsiz, -1*fft_sign, stream )
    !
    compare_len = dfft%nr2x * dfft%nr1w(me2) * dfft%my_nr3p
    IF (ny > 1) compare_len = dfft%nr2x * dfft%nr1w_tg * dfft%my_nr3p
    !
    aux(1:compare_len) = scatter_out_d(1:compare_len)
    ! Check
    CALL test%assert_close( scatter_out(1:compare_len), aux(1:compare_len) )
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
    integer :: fft_sign = 2
    integer :: vsiz, compare_len, my_nr1p_
    !
    parallel = mp%n .gt. 1
    CALL fft_desc_init(dfft, smap, "wave", gamma_only, parallel, mp%comm, nyfft=ny)
    vsiz     = dfft%nnr
    my_nr1p_ = count(dfft%ir1w > 0)
    if (ny > 1) then
       fft_sign = 3
       vsiz = dfft%nnr_tg
       my_nr1p_ = count(dfft%ir1w_tg > 0)
    end if
    !
    ! Allocate variables
    ALLOCATE(scatter_in(vsiz), scatter_out(vsiz), aux(vsiz))
    ALLOCATE(scatter_in_d(vsiz), scatter_out_d(vsiz))
    !
    ! Test 1
    CALL fill_random(scatter_in, scatter_in_d, vsiz)
    !
    CALL fft_scatter_yz( dfft, scatter_in, scatter_out, vsiz, fft_sign )
    CALL fft_scatter_yz_gpu( dfft, scatter_in_d, scatter_out_d, vsiz, fft_sign )
    ! Set the number of elements that should be strictly equivalent in the
    ! two implementations.
    compare_len = dfft%my_nr3p*my_nr1p_*dfft%nr2x
    aux(1:compare_len) = scatter_out_d(1:compare_len)
    ! Check
    CALL test%assert_close( scatter_out(1:compare_len), aux(1:compare_len) )
    !
    ! Test 2
    CALL fill_random(scatter_in, scatter_in_d, vsiz)
    !
    CALL fft_scatter_yz( dfft, scatter_out, scatter_in, vsiz, -1*fft_sign )
    CALL fft_scatter_yz_gpu( dfft, scatter_out_d, scatter_in_d, vsiz, -1*fft_sign )
    !
    compare_len = dfft%nsw(mp%me+1)*dfft%nr3x
    aux(1:compare_len) = scatter_out_d(1:compare_len)
    ! Check
    CALL test%assert_close( scatter_out(1:compare_len), aux(1:compare_len) )
    !
    CALL fft_desc_finalize(dfft, smap)
    DEALLOCATE(scatter_in, scatter_out, aux, scatter_in_d, scatter_out_d)
    !
  END SUBROUTINE test_fft_scatter_yz_gpu_1
  !
  SUBROUTINE test_fft_scatter_many_yz_gpu_1(mp, test, gamma_only, ny)
    !
    ! This test checks wave fft scatter, with parallel = .true. if 
    !  called with more than 1 MPI.
    !
    USE cudafor
    USE fft_param,       ONLY : DP
    USE fft_types,       ONLY : fft_type_descriptor
    USE stick_base,      ONLY : sticks_map
    USE scatter_mod,     ONLY : fft_scatter_yz
    USE scatter_mod_gpu, ONLY : fft_scatter_yz_gpu, fft_scatter_many_yz_gpu
    implicit none
    TYPE(mpi_t) :: mp
    TYPE(tester_t) :: test
    !
    TYPE(fft_type_descriptor) :: dfft
    TYPE(sticks_map) :: smap
    LOGICAL, INTENT(IN) :: gamma_only
    INTEGER, INTENT(IN) :: ny
    INTEGER, PARAMETER :: howmany = 4
    !
    LOGICAL :: parallel
    COMPLEX(DP), ALLOCATABLE :: scatter_in(:), scatter_in_cpy(:), scatter_out(:), aux(:)
    COMPLEX(DP), ALLOCATABLE, DEVICE :: scatter_in_d(:), scatter_out_d(:), aux_d(:)
    !                 convenient variables for slices
    integer :: i, l, start_in, end_in, start_out, end_out, nstick_zx, n3, n3x, vsiz
    integer :: start_sl, end_sl
    !integer(kind = cuda_stream_kind) :: streams(5)
    !
    parallel = mp%n .gt. 1
    IF (ny > 1) print *, 'scatter_many does not support task grouping'
    CALL fft_desc_init(dfft, smap, "wave", gamma_only, parallel, mp%comm, nyfft=ny)
    !
    ! Allocate variables
    vsiz = dfft%nnr*howmany
    ALLOCATE(scatter_in(vsiz), scatter_in_cpy(vsiz), scatter_out(vsiz), aux(vsiz))
    ALLOCATE(scatter_in_d(vsiz), scatter_out_d(vsiz), aux_d(vsiz))
    !
    ! How FFT allocates bunches in this case:
    nstick_zx = MAXVAL(dfft%nsw)
    n3 = dfft%nr3
    n3x = dfft%nr3x
    
    ! Test 1
    CALL fill_random(scatter_in, scatter_in_d, vsiz)
    scatter_in_cpy = scatter_in
    !
    !print *, 'dfft%nnr, nr3, nr2, nr1 : ', dfft%nnr, dfft%nr3, dfft%nr2, dfft%nr1 
    !  Test 1.1, compare scatter of slices 
    DO i=0,howmany-1
       start_in = i*dfft%nnr + 1
       end_in = (i+1)*dfft%nnr
       start_out = start_in
       !
       end_out = (start_out-1) + dfft%my_nr3p*dfft%nr1w(dfft%mype2 +1)*dfft%nr2x
       !
       CALL fft_scatter_yz( dfft, scatter_in(start_in:end_in), scatter_out(start_out:end_out), dfft%nnr, 2 )
       CALL fft_scatter_yz_gpu( dfft, scatter_in_d(start_in:end_in), scatter_out_d(start_in:end_out), dfft%nnr, 2 )
       aux(start_out:end_out) = scatter_out_d(start_out:end_out)
       ! Check
       CALL test%assert_close( scatter_out(start_out:end_out), aux(start_out:end_out) )
    END DO
    !
    scatter_in = scatter_in_cpy
    ! Store data as expected in input
    DO i=0,howmany-1
       start_in = i*dfft%nnr + 1
       end_in   = (i+1)*dfft%nnr
       start_out= i*nstick_zx*n3x+1
       end_out  = (i+1)*nstick_zx*n3
       scatter_in_d( start_out : end_out ) = scatter_in(start_in:start_in+nstick_zx*n3)
    END DO
    CALL fft_scatter_many_yz_gpu ( dfft, scatter_in_d, scatter_out_d, vsiz, 2, howmany )
    
    DO i=0,howmany-1
       start_out = i*dfft%nnr + 1
       end_out = (start_out-1) + dfft%my_nr3p*dfft%nr1w(dfft%mype2 +1)*dfft%nr2x
       !
       aux(start_out:end_out) = scatter_out_d(start_out:end_out)
       !
       CALL test%assert_close( scatter_out(start_out:end_out), aux(start_out:end_out) )
    END DO
    !
    !
    ! Test 2
    CALL fill_random(scatter_in, scatter_in_d, vsiz)
    scatter_in_cpy = scatter_in
    !
    DO i=0,howmany-1
       ! Input data for fft_scatter_yz call
       start_in = i*dfft%nnr + 1
       end_in = (i+1)*dfft%nnr
       !
       ! Where to store output data
       start_out = start_in
       end_out = end_in
       CALL fft_scatter_yz( dfft, scatter_out(start_out:end_out), scatter_in(start_in:end_in),  dfft%nnr, -2 )
       CALL fft_scatter_yz_gpu( dfft, scatter_out_d(start_out:end_out), scatter_in_d(start_in:end_in), dfft%nnr, -2 )
       !
       aux(start_out:end_out) = scatter_out_d(start_out:end_out)
       !
       ! Check only relevant part
       end_out = start_out + dfft%nsw(mp%me+1)*n3
       CALL test%assert_close( scatter_out(start_out:end_out), aux(start_out:end_out) )
    END DO
    !
    !
    ! Now repeat the test, but doing all the FFTs in a single shot
    scatter_in_d(1:vsiz) = scatter_in_cpy(1:vsiz)
    !
    CALL fft_scatter_many_yz_gpu ( dfft, scatter_out_d, scatter_in_d, vsiz, -2, howmany )
    
    DO i=0,howmany-1
       ! Extract data from GPU. Data are spaced by nstick_zx*n3x
       start_out = 1 + i*nstick_zx*n3x
       end_out   = (i+1)*nstick_zx*n3x
       aux(start_out:end_out) = scatter_out_d(start_out:end_out)
       !
       start_in = i*dfft%nnr + 1
       end_in   = i*dfft%nnr + dfft%nsw(mp%me+1)*n3x    !nstick_zx*nx3!
       !
       ! Extract only data tofft_scatter_yz compare the two methods. Results from the
       ! previous call to fft_scatter_yz are separated by nnr, while the
       ! new results are separated by nstick_zx*n3x. We read start_out
       ! and add from there dfft%nsw(mp%me+1)*n3x to be compared (don't forget the -1!)
       start_sl = start_out
       end_sl   = (start_out - 1) + dfft%nsw(mp%me+1)*n3x
       CALL test%assert_close( aux(start_sl:end_sl), scatter_out(start_in:end_in) )
    END DO
    !
    CALL fft_desc_finalize(dfft, smap)
    DEALLOCATE(scatter_in, scatter_in_cpy, scatter_out, aux, scatter_in_d, scatter_out_d)
    !
  END SUBROUTINE test_fft_scatter_many_yz_gpu_1
  
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
