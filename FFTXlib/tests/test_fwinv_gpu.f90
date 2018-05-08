#if defined(__CUDA)
program test_fwinv_gpu
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
    test%tolerance64 = 1.d-12
    !
    CALL save_random_seed("test_fwinv_gpu", mp%me)
    !
    DO i = 1, mp%n
      IF (MOD(mp%n,i) == 0 ) THEN
        CALL test_fwfft_gpu_1(mp, test, .true., i)
        CALL test_fwfft_gpu_1(mp, test, .false., i)
        !
        CALL test_invfft_gpu_1(mp, test, .true., i)
        CALL test_invfft_gpu_1(mp, test, .false., i)
      END IF
    END DO
    CALL test_fwfft_many_gpu_1(mp, test, .true., 1)
    CALL test_fwfft_many_gpu_1(mp, test, .false., 1)
    !
    CALL test_invfft_many_gpu_1(mp, test, .true., 1)
    CALL test_invfft_many_gpu_1(mp, test, .false., 1)
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
  !
  SUBROUTINE calc_bg(at, bg)
    USE fft_param, ONLY : DP
    implicit none
    REAL(DP), PARAMETER :: pi=4.D0*DATAN(1.D0)
    REAL(DP), INTENT(IN) :: at(3,3)
    REAL(DP), INTENT(OUT) :: bg(3,3)
    REAL(DP) :: ucvol
    !
    ucvol=at(1,1)*(at(2,2)*at(3,3)-at(3,2)*at(2,3))+&
     & at(2,1)*(at(3,2)*at(1,2)-at(1,2)*at(3,3))+&
     at(3,1)*(at(1,2)*at(2,3)-at(2,2)*at(1,3))
     
    !calculate reciprocal-space lattice vectors
    bg(1,1)=2.0*pi*(at(2,2)*at(3,3)-at(3,2)*at(2,3))/ucvol
    bg(2,1)=2.0*pi*(at(3,2)*at(1,3)-at(1,2)*at(3,3))/ucvol
    bg(3,1)=2.0*pi*(at(1,2)*at(2,3)-at(2,2)*at(1,3))/ucvol
    bg(1,2)=2.0*pi*(at(2,3)*at(3,1)-at(3,3)*at(2,1))/ucvol
    bg(2,2)=2.0*pi*(at(3,3)*at(1,1)-at(1,3)*at(3,1))/ucvol
    bg(3,2)=2.0*pi*(at(1,3)*at(2,1)-at(2,3)*at(1,1))/ucvol
    bg(1,3)=2.0*pi*(at(2,1)*at(3,2)-at(3,1)*at(2,2))/ucvol
    bg(2,3)=2.0*pi*(at(3,1)*at(1,2)-at(1,1)*at(3,2))/ucvol
    bg(3,3)=2.0*pi*(at(1,1)*at(2,2)-at(2,1)*at(1,2))/ucvol
  END SUBROUTINE calc_bg
  !
  SUBROUTINE fft_desc_init(dfft, smap, flavor, gamma_only, parallel, comm, nyfft)
    USE stick_base
    USE fft_types, ONLY : fft_type_descriptor, fft_type_init
    USE fft_param, ONLY : DP
    implicit none
    TYPE(fft_type_descriptor) :: dfft
    TYPE(sticks_map) :: smap
    CHARACTER(LEN=*), INTENT(IN) :: flavor
    LOGICAL :: gamma_only
    LOGICAL :: parallel
    INTEGER :: comm, nyfft
    !
    REAL(DP) :: at(3,3), bg(3,3)
    !
    at = RESHAPE((/10.d0, 0.d0, 0.d0, 0.d0, 10.d0, 0.d0, 0.d0, 0.d0, 10.d0/), shape(at))
    CALL calc_bg(at, bg)
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
  SUBROUTINE test_fwfft_gpu_1(mp, test, gamma_only, ny)
    USE cudafor
    USE fft_param,       ONLY : DP
    USE fft_types,       ONLY : fft_type_descriptor
    USE stick_base,      ONLY : sticks_map
    USE fft_interfaces,  ONLY : fwfft
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
    COMPLEX(DP), ALLOCATABLE :: data_in(:), aux(:)
    COMPLEX(DP), ALLOCATABLE, DEVICE :: data_in_d(:)
    INTEGER :: i
    !
    parallel = mp%n .gt. 1
    CALL fft_desc_init(dfft, smap, 'wave', gamma_only, parallel, mp%comm, nyfft=ny)
    dfft%rho_clock_label='bla' ; dfft%wave_clock_label='bla'
    !
    ! Test 1
    !
    IF ( ny .gt. 1 ) THEN
      ! Allocate variables
      ALLOCATE(data_in(dfft%nnr_tg), aux(dfft%nnr_tg))
      ALLOCATE(data_in_d(dfft%nnr_tg))
      CALL fill_random(data_in, data_in_d, dfft%nnr_tg)
      !
      CALL fwfft( 'tgWave' , data_in, dfft, 1 )
      CALL fwfft( 'tgWave' , data_in_d, dfft, 1 )
    ELSE
      ALLOCATE(data_in(dfft%nnr), aux(dfft%nnr))
      ALLOCATE(data_in_d(dfft%nnr))
      CALL fill_random(data_in, data_in_d, dfft%nnr)
      !
      CALL fwfft( 'Wave' , data_in, dfft, 1 )
      CALL fwfft( 'Wave' , data_in_d, dfft, 1 )    
    ENDIF
    aux = data_in_d
    ! Check
    CALL test%assert_close( data_in(1:dfft%ngw), aux(1:dfft%ngw) )
    !
    ! Test 2
    !
    DEALLOCATE(data_in, data_in_d, aux)
    ALLOCATE(data_in(dfft%nnr), aux(dfft%nnr))
    ALLOCATE(data_in_d(dfft%nnr))
    CALL fill_random(data_in, data_in_d, dfft%nnr)
    !
    CALL fwfft( 'Rho' , data_in, dfft, 1 )
    CALL fwfft( 'Rho' , data_in_d, dfft, 1 )
    aux = data_in_d
    ! Check
    CALL test%assert_close( data_in, aux )
    !
    CALL fft_desc_finalize(dfft, smap)
    DEALLOCATE(data_in, data_in_d, aux)
    !
  END SUBROUTINE test_fwfft_gpu_1
  !
  SUBROUTINE test_invfft_gpu_1(mp, test, gamma_only, ny)
    USE cudafor
    USE fft_param,       ONLY : DP
    USE fft_types,       ONLY : fft_type_descriptor
    USE stick_base,      ONLY : sticks_map
    USE fft_interfaces,  ONLY : invfft
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
    COMPLEX(DP), ALLOCATABLE :: data_in(:), aux(:)
    COMPLEX(DP), ALLOCATABLE, DEVICE :: data_in_d(:)
    !
    parallel = mp%n .gt. 1
    CALL fft_desc_init(dfft, smap, 'wave', gamma_only, parallel, mp%comm, nyfft=ny)
    dfft%rho_clock_label='bla' ; dfft%wave_clock_label='bla'
    !
    ! Test 1
    !
    IF ( ny .gt. 1 ) THEN
      ! Allocate variables
      ALLOCATE(data_in(dfft%nnr_tg), aux(dfft%nnr_tg))
      ALLOCATE(data_in_d(dfft%nnr_tg))
      CALL fill_random(data_in, data_in_d, dfft%nnr_tg)
      !
      CALL invfft( 'tgWave' , data_in, dfft, 1 )
      CALL invfft( 'tgWave' , data_in_d, dfft, 1 )
    ELSE
      !
      ! Allocate variables
      ALLOCATE(data_in(dfft%nnr), aux(dfft%nnr))
      ALLOCATE(data_in_d(dfft%nnr))
      CALL fill_random(data_in, data_in_d, dfft%nnr)
      !
      CALL invfft( 'Wave' , data_in, dfft, 1 )
      CALL invfft( 'Wave' , data_in_d, dfft, 1 )    
    ENDIF
    aux = data_in_d
    ! Check
    CALL test%assert_close( data_in, aux )
    !
    ! Test 2
    !
    DEALLOCATE(data_in, data_in_d, aux)
    ALLOCATE(data_in(dfft%nnr), aux(dfft%nnr))
    ALLOCATE(data_in_d(dfft%nnr))
    CALL fill_random(data_in, data_in_d, dfft%nnr)
    !
    CALL invfft( 'Rho' , data_in, dfft, 1 )
    CALL invfft( 'Rho' , data_in_d, dfft, 1 )
    aux = data_in_d
    ! Check
    CALL test%assert_close( data_in, aux )
    !
    CALL fft_desc_finalize(dfft, smap)
    DEALLOCATE(data_in, data_in_d, aux)
    !
  END SUBROUTINE test_invfft_gpu_1
  !
  SUBROUTINE test_fwfft_many_gpu_1(mp, test, gamma_only, ny)
    USE cudafor
    USE fft_param,       ONLY : DP
    USE fft_types,       ONLY : fft_type_descriptor
    USE stick_base,      ONLY : sticks_map
    USE fft_interfaces,  ONLY : fwfft
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
    COMPLEX(DP), ALLOCATABLE :: data_in(:), aux(:)
    COMPLEX(DP), ALLOCATABLE, DEVICE :: data_in_d(:)
    integer, parameter :: howmany=4
    INTEGER :: i, start
    !
    parallel = mp%n .gt. 1
    CALL fft_desc_init(dfft, smap, 'wave', gamma_only, parallel, mp%comm, nyfft=ny)
    dfft%rho_clock_label='bla' ; dfft%wave_clock_label='bla'
    !
    ! Test 1
    !
    IF ( ny .gt. 1 ) THEN
      ! Not (yet?) possible
      RETURN
    ELSE
      ALLOCATE(data_in(dfft%nnr*howmany), aux(dfft%nnr*howmany))
      ALLOCATE(data_in_d(dfft%nnr*howmany))
      CALL fill_random(data_in, data_in_d, dfft%nnr*howmany)
      !
      CALL fwfft( 'Wave' , data_in_d, dfft, howmany=howmany)
      !
      DO i=0,0
        start = i*dfft%nnr
        CALL fwfft( 'Wave' , data_in(1+start:), dfft, 1 )
        aux(start+1:start+dfft%nnr) = data_in_d(start+1:start+dfft%nnr)
        ! Check
        CALL test%assert_close( data_in(start+1:start+dfft%ngw), aux(start+1:start+dfft%ngw) )
      END DO
      !
    ENDIF
    !
    ! Test 2
    !
    !!!! DEALLOCATE(data_in, data_in_d, aux)
    !!!! ALLOCATE(data_in(dfft%nnr), aux(dfft%nnr))
    !!!! ALLOCATE(data_in_d(dfft%nnr))
    !!!! CALL fill_random(data_in, data_in_d, dfft%nnr)
    !!!! !
    !!!! CALL fwfft( 'Rho' , data_in, dfft, 1 )
    !!!! CALL fwfft( 'Rho' , data_in_d, dfft, 1 )
    !!!! aux = data_in_d
    !!!! ! Check
    !!!! CALL test%assert_close( data_in, aux )
    !
    CALL fft_desc_finalize(dfft, smap)
    DEALLOCATE(data_in, data_in_d, aux)
    !
  END SUBROUTINE test_fwfft_many_gpu_1
  !
  SUBROUTINE test_invfft_many_gpu_1(mp, test, gamma_only, ny)
    USE cudafor
    USE fft_param,       ONLY : DP
    USE fft_types,       ONLY : fft_type_descriptor
    USE stick_base,      ONLY : sticks_map
    USE fft_interfaces,  ONLY : invfft
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
    COMPLEX(DP), ALLOCATABLE :: data_in(:), aux(:)
    COMPLEX(DP), ALLOCATABLE, DEVICE :: data_in_d(:)
    INTEGER(kind = cuda_stream_kind) :: strm = 0
    integer, parameter :: howmany=4
    integer :: start, i
    !
    parallel = mp%n .gt. 1
    CALL fft_desc_init(dfft, smap, 'wave', gamma_only, parallel, mp%comm, nyfft=ny)
    dfft%rho_clock_label='bla' ; dfft%wave_clock_label='bla'
    !
    ! Test 1
    !
    IF ( ny .gt. 1 ) THEN
      ! Not (yet?) possible
      RETURN
    ELSE
      !
      ! Allocate variables
      ALLOCATE(data_in(howmany*dfft%nnr), aux(howmany*dfft%nnr))
      ALLOCATE(data_in_d(howmany*dfft%nnr))
      CALL fill_random(data_in, data_in_d, howmany*dfft%nnr)
      !
      !CALL invfft( 'Wave' , data_in, dfft, 1 )
      CALL invfft( 'Wave' , data_in_d, dfft, howmany=howmany, stream=strm )
      DO i=0,howmany-1
        start = i*dfft%nnr
        CALL invfft( 'Wave' , data_in(1+start:), dfft, 1 )
        aux(start+1:start+dfft%nnr) = data_in_d(start+1:start+dfft%nnr)
        ! Check
        CALL test%assert_close( data_in(start+1:start+dfft%nnr), aux(start+1:start+dfft%nnr) )
      END DO
    ENDIF
    !aux = data_in_d
    ! Check
    !CALL test%assert_close( data_in, aux )
    !
    ! Test 2
    !
    !!!!  DEALLOCATE(data_in, data_in_d, aux)
    !!!!  ALLOCATE(data_in(dfft%nnr), aux(dfft%nnr))
    !!!!  ALLOCATE(data_in_d(dfft%nnr))
    !!!!  CALL fill_random(data_in, data_in_d, dfft%nnr)
    !!!!  !
    !!!!  CALL invfft( 'Rho' , data_in, dfft, 1 )
    !!!!  CALL invfft( 'Rho' , data_in_d, dfft, 1 )
    !!!!  aux = data_in_d
    !!!!  ! Check
    !!!!  CALL test%assert_close( data_in, aux )
    !!!!  !
    !!!!  CALL fft_desc_finalize(dfft, smap)
    !!!!  DEALLOCATE(data_in, data_in_d, aux)
    !
  END SUBROUTINE test_invfft_many_gpu_1
  
end program test_fwinv_gpu
! dummy subroutines
subroutine stop_clock( label )
character(len=*) :: label
end subroutine stop_clock
subroutine start_clock( label )
character(len=*) :: label
end subroutine start_clock
!
#else
program test_fwinv_gpu
end program test_fwinv_gpu
#endif
