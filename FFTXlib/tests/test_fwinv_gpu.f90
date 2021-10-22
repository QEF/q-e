#if defined(__CUDA)
#define DEVATTR ,DEVICE
#else
#define DEVATTR
#endif

program test_fwinv_gpu
#if defined(__MPI) && defined(__MPI_MODULE)
    USE mpi
#endif
    USE tester
    IMPLICIT NONE
#if defined(__MPI) && ! defined(__MPI_MODULE)
    INCLUDE 'mpif.h'
#endif
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
    ! Prepare MPI and communicators
    CALL mpi_data_init(mp%me, mp%n, mp%root, mp%comm)
    !
    CALL test%init()
    !
    ! A rather large threshold is necessary to match the results of all
    ! possible implementations.
    test%tolerance64 = 1.d-10
    !
    CALL save_random_seed("test_fwinv_gpu", mp%me)
    !
    DO i = 1, mp%n
      IF (MOD(mp%n,i) == 0 ) THEN
        ! test R -> G, gamma case
        CALL test_fwfft_gpu_1(mp, test, .true., i)
        ! test R -> G, k case
        CALL test_fwfft_gpu_1(mp, test, .false., i)
        !
        ! test G -> R, gamma case
        CALL test_invfft_gpu_1(mp, test, .true., i)
        ! test G -> R, k case
        CALL test_invfft_gpu_1(mp, test, .false., i)
      END IF
    END DO
    !
#if defined(__CUDA)
    ! the batched FFT is only implemented for GPU,
    ! will segault if called on non-implemented version.
    CALL test_fwfft_many_gpu_1(mp, test, .true., 1)
    CALL test_fwfft_many_gpu_1(mp, test, .false., 1)
    !
    CALL test_invfft_many_gpu_1(mp, test, .true., 1)
    CALL test_invfft_many_gpu_1(mp, test, .false., 1)
#endif
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
    INTEGER :: ngm, ngm_g
    !
    REAL(DP) :: at(3,3), bg(3,3), alat, tpiba
    REAL(DP), ALLOCATABLE :: g(:, :)
    !
    REAL(DP), PARAMETER :: gcut = 80.d0
    REAL(DP), PARAMETER :: dual = 4.d0
    REAL(DP), PARAMETER :: pi=4.D0*DATAN(1.D0)
    !
    ! Define a direct lattice
    at = RESHAPE((/10.d0, 0.d0, 0.d0, 0.d0, 10.d0, 0.d0, 0.d0, 0.d0, 10.d0/), shape(at))
    !
    alat = SQRT ( at(1,1)**2+at(2,1)**2+at(3,1)**2 )
    !
    ! Lattice must be defined in units of alat
    at(:,:) = at(:,:) / alat
    !
    tpiba = 2.0d0*pi/alat
    !
    ! And the related recuprocal space
    CALL recips(at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
    !
    ! In a FFT of flavor='wave' the dual, here set to 4.0d0, will multiply gcut to obtain
    ! the cutoff for hosting "charges" i.e. wfc**2.
    !
    ! Inputs are:
    ! dfft       -> the fft type containing all details about the distributed grid
    ! smap       -> the map of the sticks used to distribute the data in rec. space
    ! flavor     -> can be be 'wave' or 'rho'
    ! gamma_only -> setup the complex FFT to perform two real FFTs at a time
    ! parallel   -> true if more than 1 process will take part in the FFT, false otherwise
    ! comm       -> the mpi communicator to be used for messages among the processes performing the FFTs,
    !                basically alltoall communications to perform (partial) transposition of the data
    ! at         -> the real space lattice, needs to be define in units of alat
    ! bg         -> the reciprocal lattice
    ! gcut       -> the cutoff for the plane wave expansion for flavor='wave'
    ! dual       -> the dimension of the sphere of g vectors (generally meant
    !                to be large enough to host products of wavefunctions)
    ! nyfft      -> data are decomposed along Z in real space when parallel=.true.,
    !                along X in real space when parallel=.false., and as sticks in reciprocal
    !                space. nyfft describes a second level of division along another dimension in real space
    !                or, when enabled, as portions of multiple bands. ???
    ! nmany      -> maximum number of bands to be transformed simultaneously ???
    !
    CALL fft_type_init(dfft, smap, flavor, gamma_only, parallel, comm, at, bg, gcut, dual, &
    & nyfft=nyfft, nmany=1)
    !
    ! Create a map between g vectors and distributed grid point in the fft.
    ! We will not he gvectors, but the map, stored in dfft$nl and dfft%nlm (nl minus, for gamma case)
    ! will be used to check only physically relevant numbers.
    !
    ALLOCATE (g(3, dfft%ngm))
    !
    ! Largest g vector size
    ngm = dfft%ngm
    !
#if defined(__MPI)
    CALL MPI_ALLREDUCE(ngm, ngm_g, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
    ngm_g = ngm
#endif
    
    ! Generate G vectors and map global g vectors to local FFT points.
    CALL ggen(dfft, gamma_only, at, bg, 4.d0*gcut, ngm_g, ngm, g, .false.)
    !
    DEALLOCATE(g)
    !
  END SUBROUTINE fft_desc_init
  !
  SUBROUTINE ggen ( dfftp, gamma_only, at, bg,  gcutm, ngm_g, ngm, g, no_global_sort )
    !----------------------------------------------------------------------
    !
    !     This routine generates all the reciprocal lattice vectors
    !     contained in the sphere of radius gcutm. Furthermore it
    !     computes the indices nl which give the correspondence
    !     between the fft mesh points and the array of g vectors.
    !
    USE fft_types, ONLY: fft_stick_index, fft_type_descriptor
    USE fft_ggen, ONLY : fft_set_nl
    USE fft_param
    !USE mp, ONLY: mp_rank, mp_size, mp_sum
    !USE constants, ONLY : eps8
    !
    IMPLICIT NONE
    !
    TYPE(fft_type_descriptor),INTENT(INOUT) :: dfftp
    LOGICAL,  INTENT(IN) :: gamma_only
    REAL(DP), INTENT(IN) :: at(3,3), bg(3,3), gcutm
    INTEGER, INTENT(IN) :: ngm_g
    INTEGER, INTENT(INOUT) :: ngm
    REAL(DP), INTENT(OUT) :: g(:,:)
    !  if no_global_sort is present (and it is true) G vectors are sorted only
    !  locally and not globally. In this case no global array needs to be
    !  allocated and sorted: saves memory and a lot of time for large systems.
    !
    LOGICAL, INTENT(IN) :: no_global_sort
    !
    !     here a few local variables
    !
    REAL(DP) :: tx(3), ty(3), t(3)
    REAL(DP), ALLOCATABLE :: tt(:)
    INTEGER :: ngm_save, ngm_max, ngm_local
    !
    REAL(DP), ALLOCATABLE :: g2sort_g(:)
    ! array containing only g vectors for the current processor
    INTEGER, ALLOCATABLE :: mill_unsorted(:,:)
    ! array containing all g vectors generators, on all processors
    ! (replicated data). When no_global_sort is present and .true.,
    ! only g-vectors for the current processor are stored
    INTEGER, ALLOCATABLE :: igsrt(:), g2l(:)
    !
    INTEGER :: ni, nj, nk, i, j, k, ng
    INTEGER :: istart, jstart, kstart
    LOGICAL :: global_sort, is_local
    !
    ! The 'no_global_sort' is not optional in this case.
    ! This differs from the version present in QE distribution.
    global_sort = .NOT. no_global_sort
    !
    IF( .NOT. global_sort ) THEN
       ngm_max = ngm
    ELSE
       ngm_max = ngm_g
    END IF
    !
    ! save current value of ngm
    !
    ngm_save  = ngm
    !
    ngm = 0
    ngm_local = 0
    !
    !    set the total number of fft mesh points and initial value of gg
    !    The choice of gcutm is due to the fact that we have to order the
    !    vectors after computing them.
    !
    !
    !    and computes all the g vectors inside a sphere
    !
    ALLOCATE( mill_unsorted( 3, ngm_save ) )
    ALLOCATE( igsrt( ngm_max ) )
    ALLOCATE( g2l( ngm_max ) )
    ALLOCATE( g2sort_g( ngm_max ) )
    !
    g2sort_g(:) = 1.0d20
    !
    ! allocate temporal array
    !
    ALLOCATE( tt( dfftp%nr3 ) )
    !
    ! max miller indices (same convention as in module stick_set)
    !
    ni = (dfftp%nr1-1)/2
    nj = (dfftp%nr2-1)/2
    nk = (dfftp%nr3-1)/2
    !
    ! gamma-only: exclude space with x < 0
    !
    IF ( gamma_only ) THEN
       istart = 0
    ELSE
       istart = -ni
    ENDIF
    !
    iloop: DO i = istart, ni
       !
       ! gamma-only: exclude plane with x = 0, y < 0
       !
       IF ( gamma_only .and. i == 0 ) THEN
          jstart = 0
       ELSE
          jstart = -nj
       ENDIF
       !
       tx(1:3) = i * bg(1:3,1)
       !
       jloop: DO j = jstart, nj
          !
          IF ( .NOT. global_sort ) THEN
             IF ( fft_stick_index( dfftp, i, j ) == 0 ) CYCLE jloop
             is_local = .TRUE.
          ELSE
             IF ( dfftp%lpara .AND. fft_stick_index( dfftp, i, j ) == 0) THEN
                is_local = .FALSE.
             ELSE
                is_local = .TRUE.
             END IF
          END IF
          !
          ! gamma-only: exclude line with x = 0, y = 0, z < 0
          !
          IF ( gamma_only .and. i == 0 .and. j == 0 ) THEN
             kstart = 0
          ELSE
             kstart = -nk
          ENDIF
          !
          ty(1:3) = tx(1:3) + j * bg(1:3,2)
          !
          !  compute all the norm square
          !
          DO k = kstart, nk
             !
             t(1) = ty(1) + k * bg(1,3)
             t(2) = ty(2) + k * bg(2,3)
             t(3) = ty(3) + k * bg(3,3)
             tt(k-kstart+1) = t(1)**2 + t(2)**2 + t(3)**2
          ENDDO
          !
          !  save all the norm square within cutoff
          !
          DO k = kstart, nk
             IF (tt(k-kstart+1) <= gcutm) THEN
                ngm = ngm + 1
                IF (ngm > ngm_max) CALL fftx_error__ ('ggen 1', 'too many g-vectors', ngm)
                IF ( tt(k-kstart+1) > eps8 ) THEN
                   g2sort_g(ngm) = tt(k-kstart+1)
                ELSE
                   g2sort_g(ngm) = 0.d0
                ENDIF
                IF (is_local) THEN
                  ngm_local = ngm_local + 1
                  mill_unsorted( :, ngm_local ) = (/ i,j,k /)
                  g2l(ngm) = ngm_local
                ELSE
                  g2l(ngm) = 0
                ENDIF
             ENDIF
          ENDDO
       ENDDO jloop
    ENDDO iloop

    IF (ngm  /= ngm_max) &
         CALL fftx_error__ ('ggen', 'g-vectors missing !', abs(ngm - ngm_max))
    !
    igsrt(1) = 0
    IF( .NOT. global_sort ) THEN
       CALL hpsort_eps( ngm, g2sort_g, igsrt, eps8 )
    ELSE
       CALL hpsort_eps( ngm_g, g2sort_g, igsrt, eps8 )
    END IF
    DEALLOCATE( g2sort_g, tt )
    !
    ngm = 0
    !
    ngloop: DO ng = 1, ngm_max
       !
       IF (g2l(igsrt(ng))>0) THEN
          ! fetch the indices
          i = mill_unsorted(1, g2l(igsrt(ng)))
          j = mill_unsorted(2, g2l(igsrt(ng)))
          k = mill_unsorted(3, g2l(igsrt(ng)))
          !
          ngm = ngm + 1
          !
          !  Here map local and global g index !!! N.B: :
          !  the global G vectors arrangement depends on the number of processors
          !
          g(1:3, ngm) = i * bg (:, 1) + j * bg (:, 2) + k * bg (:, 3)
          !
       ENDIF
    ENDDO ngloop

    DEALLOCATE( igsrt, g2l )
    IF (ngm /= ngm_save) &
         CALL fftx_error__ ('ggen', 'g-vectors (ngm) missing !', abs(ngm - ngm_save))
    !
    !     Now set nl and nls with the correct fft correspondence
    !
    CALL fft_set_nl( dfftp, at, g)
    !
  END SUBROUTINE ggen
  !
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
    USE fft_param, ONLY : DP
    implicit none
    complex(DP) DEVATTR :: c_d(:)
    complex(DP)         :: c(:)
    integer, intent(in) :: n
    !
    real(DP), ALLOCATABLE :: rnd_aux(:)
    !
    ALLOCATE (rnd_aux(2*n))
    CALL RANDOM_NUMBER(rnd_aux)
    c(1:n) = CMPLX(rnd_aux(1:n), rnd_aux(n+1:2*n), kind=DP)
    c_d = c
    DEALLOCATE(rnd_aux)
  END SUBROUTINE fill_random
  !
  SUBROUTINE fill_random_cpu(c, n)
    USE fft_param, ONLY : DP
    implicit none
    complex(DP)         :: c(:)
    integer, intent(in) :: n
    !
    real(DP), ALLOCATABLE :: rnd_aux(:)
    !
    ALLOCATE (rnd_aux(2*n))
    CALL RANDOM_NUMBER(rnd_aux)
    c(1:n) = CMPLX(rnd_aux(1:n), rnd_aux(n+1:2*n), kind=DP)
    DEALLOCATE(rnd_aux)
  END SUBROUTINE fill_random_cpu
  !
  SUBROUTINE test_fwfft_gpu_1(mp, test, gamma_only, ny)
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
    COMPLEX(DP), ALLOCATABLE DEVATTR :: data_in_d(:)
    INTEGER :: i
    !
    ! task groups not implemented in 2D decomposition. Need to check the other case
    IF ( ny .gt. 1 ) return

    parallel = mp%n .gt. 1
    CALL fft_desc_init(dfft, smap, 'wave', gamma_only, parallel, mp%comm, nyfft=ny)
    dfft%rho_clock_label='bla' ; dfft%wave_clock_label='bla'
    !
    ! Test 1
    !
    IF ( ny .gt. 1 ) THEN
      ! Allocate variables and fill realspace data with random numbers
      ALLOCATE(data_in(dfft%nnr_tg), aux(dfft%nnr_tg))
      ALLOCATE(data_in_d(dfft%nnr_tg))
      CALL fill_random(data_in, data_in_d, dfft%nnr_tg)
      !
      CALL fwfft( 'tgWave' , data_in, dfft, 1 )
      CALL fwfft( 'tgWave' , data_in_d, dfft, 1 )
    ELSE
      ! Allocate variables and fill realspace data with random numbers
      ALLOCATE(data_in(dfft%nnr), aux(dfft%nnr))
      ALLOCATE(data_in_d(dfft%nnr))
      CALL fill_random(data_in, data_in_d, dfft%nnr) 
      !
      CALL fwfft( 'Wave' , data_in, dfft, 1 )
      CALL fwfft( 'Wave' , data_in_d, dfft, 1 )    
    ENDIF
    ! data from GPU is moved to an auxiliary array to compare the results of the GPU
    ! and the CPU implementation on the host
    aux = data_in_d
    !
    ! Check, only the values relevant for a wavefunction FFT are considered
    DO i=1,dfft%ngw
      IF (gamma_only) CALL test%assert_close( data_in(dfft%nlm(i)), aux(dfft%nlm(i)) )
      IF (.not. gamma_only) CALL test%assert_close( data_in(dfft%nl(i)), aux(dfft%nl(i)) )
    ENDDO
    !
    DEALLOCATE(data_in, data_in_d, aux)
    !
    !
    ! Test 2
    !
    ! Same as above
    ALLOCATE(data_in(dfft%nnr), aux(dfft%nnr))
    ALLOCATE(data_in_d(dfft%nnr))
    CALL fill_random(data_in, data_in_d, dfft%nnr)
    !
    CALL fwfft( 'Rho' , data_in, dfft, 1 )
    CALL fwfft( 'Rho' , data_in_d, dfft, 1 )
    aux = data_in_d
    !
    ! Check, only the values relevant for a product of wavefunctions are considered
    !
    DO i=1,dfft%ngm
      IF (gamma_only) CALL test%assert_close( data_in(dfft%nlm(i)), aux(dfft%nlm(i)) )
      IF (.not. gamma_only) CALL test%assert_close( data_in(dfft%nl(i)), aux(dfft%nl(i)) )
    ENDDO
    !
    CALL fft_desc_finalize(dfft, smap)
    DEALLOCATE(data_in, data_in_d, aux)
    !
  END SUBROUTINE test_fwfft_gpu_1
  !
  SUBROUTINE test_invfft_gpu_1(mp, test, gamma_only, ny)
    !
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
    COMPLEX(DP), ALLOCATABLE DEVATTR :: data_in_d(:)
    integer :: i
    !

    ! task groups not implemented in 2D decomposition. Need to check the other case
    IF ( ny .gt. 1 ) return

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
      !
      ! Data here is not correctly filled, but this test is disabled.
      ! This is left as TODO!!!
      !
      CALL fill_random(data_in, data_in_d, dfft%nnr_tg)
      !
      CALL invfft( 'tgWave' , data_in, dfft, 1 )
      CALL invfft( 'tgWave' , data_in_d, dfft, 1 )
    ELSE
      !
      ALLOCATE(data_in(dfft%nnr), aux(dfft%nnr))
      ALLOCATE(data_in_d(dfft%nnr))
      !
      ! Prepare input data, only vectors of wavefunctions
      data_in = (0.d0, 0.d0)
      CALL fill_random_cpu(aux, dfft%ngw)
      DO i=1, dfft%ngw
        IF (gamma_only)       data_in(dfft%nlm(i)) = aux(i)
        IF (.not. gamma_only) data_in(dfft%nl(i)) = aux(i)
      ENDDO
      ! copy to gpu and cleanup aux
      data_in_d = data_in
      aux = (0.d0, 0.d0)
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
    ! Prepare input data
    data_in = (0.d0, 0.d0)
    CALL fill_random_cpu(aux, dfft%ngm)
    DO i=1, dfft%ngm
      IF (gamma_only)       data_in(dfft%nlm(i)) = aux(i)
      IF (.not. gamma_only) data_in(dfft%nl(i)) = aux(i)
    ENDDO
    ! copy to gpu and cleanup aux
    data_in_d = data_in
    aux = (0.d0, 0.d0)
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
    !
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
    COMPLEX(DP), ALLOCATABLE DEVATTR :: data_in_d(:)
    integer, parameter :: howmany=4
    INTEGER :: i, ii, start
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
      DO i=0,howmany-1
        start = i*dfft%nnr
        CALL fwfft( 'Wave' , data_in(1+start:), dfft, 1 )
        aux(1:dfft%nnr) = data_in_d(start+1:start+dfft%nnr)
        ! Check, only the values relevant for a wavefunction FFT are considered
        DO ii=1,dfft%ngw
            IF (gamma_only) CALL test%assert_close( data_in(dfft%nlm(ii)+start), aux(dfft%nlm(ii)) )
            IF (.not. gamma_only) CALL test%assert_close( data_in(dfft%nl(ii)+start), aux(dfft%nl(ii)) )
        ENDDO
        !
      END DO
      !
    ENDIF
    !
    ! Test 2
    !
    DEALLOCATE(data_in, data_in_d, aux)
    ALLOCATE(data_in(dfft%nnr*howmany), aux(dfft%nnr))
    ALLOCATE(data_in_d(dfft%nnr*howmany))
    !
    CALL fill_random(data_in, data_in_d, dfft%nnr*howmany)
    !
    CALL fwfft( 'Rho' , data_in_d, dfft,  howmany)
    DO i=0,howmany-1
      !
      start = i*dfft%nnr
      !
      ! This will FFT the content of data_in starting from start+1 and for nnr elements
      CALL fwfft( 'Rho' , data_in(1+start:), dfft, 1 )
      aux(1:dfft%nnr) = data_in_d(start+1:start+dfft%nnr)
      !
      ! Same check as above, but remember that data_in starts from "start"
      DO ii=1,dfft%ngm
        IF (gamma_only) CALL test%assert_close( data_in(dfft%nlm(ii)+start), aux(dfft%nlm(ii)) )
        IF (.not. gamma_only) CALL test%assert_close( data_in(dfft%nl(ii)+start), aux(dfft%nl(ii)) )
      ENDDO
      !
    END DO
    !
    CALL fft_desc_finalize(dfft, smap)
    DEALLOCATE(data_in, data_in_d, aux)
    !
  END SUBROUTINE test_fwfft_many_gpu_1
  !
  SUBROUTINE test_invfft_many_gpu_1(mp, test, gamma_only, ny)
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
    COMPLEX(DP), ALLOCATABLE DEVATTR :: data_in_d(:)
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
      !
      data_in = (0.d0, 0.d0)
      CALL fill_random_cpu(aux, dfft%ngw)
      DO i=1, dfft%ngw
        IF (gamma_only)       data_in(dfft%nlm(i)) = aux(i)
        IF (.not. gamma_only) data_in(dfft%nl(i)) = aux(i)
      ENDDO
      ! copy data to simulate multiple bands
      DO i=0,howmany-1
        start = i*dfft%nnr
        data_in(start+1:start+dfft%nnr) = data_in(1:dfft%nnr)
      ENDDO
      ! copy to gpu and cleanup aux
      data_in_d = data_in
      aux = (0.d0, 0.d0)
      !
      CALL invfft( 'Wave' , data_in_d, dfft, howmany=howmany ) !, stream=strm )
      DO i=0,howmany-1
        start = i*dfft%nnr
        CALL invfft( 'Wave' , data_in(1+start:), dfft, 1 )
        aux(start+1:start+dfft%nnr) = data_in_d(start+1:start+dfft%nnr)
        ! Check
        CALL test%assert_close( data_in(start+1:start+dfft%nnr), aux(start+1:start+dfft%nnr) )
      END DO
    ENDIF
    !
    ! Test 2
    !
    DEALLOCATE(data_in, data_in_d, aux)
    ALLOCATE(data_in(dfft%nnr*howmany), aux(dfft%nnr))
    ALLOCATE(data_in_d(dfft%nnr*howmany))
    data_in = (0.d0, 0.d0)
    CALL fill_random_cpu(aux, dfft%ngm)
    !
    ! Prepare vectors assuming that a product of wfcs in reciprocal space
    ! is stored in data_in
    !
    DO i=1, dfft%ngm
      IF (gamma_only)       data_in(dfft%nlm(i)) = aux(i)
      IF (.not. gamma_only) data_in(dfft%nl(i)) = aux(i)
    ENDDO
    ! copy data to simulate multiple bands
    DO i=0,howmany-1
        start = i*dfft%nnr
        data_in(start+1:start+dfft%nnr) = data_in(1:dfft%nnr)
    ENDDO
    ! copy to gpu input data and cleanup aux
    data_in_d = data_in
    aux = (0.d0, 0.d0)

    !
    CALL invfft( 'Rho' , data_in_d, dfft, howmany )
    !
    DO i=0,howmany-1
      start = i*dfft%nnr
      CALL invfft( 'Rho' , data_in(1+start:), dfft, 1 )
      aux(1:dfft%nnr) = data_in_d(start+1:start+dfft%nnr)
      ! Check
      CALL test%assert_close( data_in(start+1:start+dfft%nnr), aux(1:dfft%nnr) )
    END DO
    !
    CALL fft_desc_finalize(dfft, smap)
    DEALLOCATE(data_in, data_in_d, aux)
    !
  END SUBROUTINE test_invfft_many_gpu_1
  
end program test_fwinv_gpu
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
