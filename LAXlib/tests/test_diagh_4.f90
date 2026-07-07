! This test uses the internal parallel diagonalization algorithm of LAXlib
! to solve standard eigenvalue problems (no overlap matrix).
!
! If the scalapack or ELPA driver is used, the test is skipped.
!
program test_diagh_4

    USE laxlib_parallel_include
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : mp_world_start, mp_world_end, mpime, &
                              root, world_comm
    USE mp_bands_util, ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
    USE tester
    USE test_helpers,  ONLY : hermitian, symmetric
    IMPLICIT NONE
    include 'laxlib_kinds.fh'
    include 'laxlib_param.fh'
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    !
    CALL test%init()
    test%tolerance64=1.d-8
    !
#if defined(__MPI)
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    !
    me_bgrp = mpime; root_bgrp=root; intra_bgrp_comm=world_comm
    !
    CALL parallel_real_1(test)
    !
    CALL parallel_complex_1(test)
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
  CONTAINS
  !
  SUBROUTINE parallel_real_1(test)
    USE mp_world,    ONLY : mpime
    USE LAXlib
    USE laxlib_descriptor, ONLY : la_descriptor, descla_init, laxlib_desc_to_intarray
    implicit none
    !
    TYPE(tester_t) :: test
    !
    TYPE(la_descriptor) :: desc
    INTEGER :: idesc(LAX_DESC_SIZE)
    integer :: n, m
    real(DP), allocatable :: h(:,:)            !< full H matrix
    real(DP), allocatable :: h_save(:,:)       !< full H, used to check consistence across calls
    real(DP), allocatable    :: e(:)           !< full set of eigenvalues
    real(DP), allocatable :: v(:,:)            !< full eigenvectors
    real(DP), allocatable    :: e_save(:)      !< full set of eigenvalues, used for checks
    real(DP), allocatable :: v_save(:,:)       !< full set of eigenvectors, used for checks
    !
    integer                  :: i, j, ii, jj, info
    logical                  :: la_proc
    !
    ! Test with a simple symmetric matrix
    n = 100
    m = n
    !
    ALLOCATE(h(n,n), e(n), v(n,n))
    ALLOCATE(h_save(n,n), e_save(n), v_save(n,n))
    !
    ! Create a random symmetric matrix on rank 0 and broadcast to all
    IF (mpime == 0) CALL symmetric(n, h)
    CALL mp_bcast(h, 0, world_comm)
    h_save = h
    !
    ! Solve with serial version first to get reference
    CALL diagh(n, m, h_save, e_save, v_save, me_bgrp, root_bgrp, intra_bgrp_comm)
    !
    ! Now test parallel version
    CALL init_parallel_diag(desc, n)
    !
    CALL laxlib_desc_to_intarray( idesc, desc )
    !
    ! Poison v so an unbroadcast rank is detectable.
    v = 1234.5_DP
    e = 0.d0
    CALL pdiagh( n, h, e, v, idesc )
    !
    CALL test%assert_close( e(1:m), e_save(1:m) )
    !
    ! Run on every rank; |<v_i, v_save_i>| = 1 catches v left
    ! uninitialized on ranks outside the ortho pool.
    DO i = 1, m
       CALL test%assert_close( ABS(SUM(v(:,i) * v_save(:,i))), 1.0_DP )
    END DO
    !
    DEALLOCATE(h, e, v, h_save, e_save, v_save)
    !
  END SUBROUTINE parallel_real_1
  !
  SUBROUTINE parallel_complex_1(test)
    USE mp_world, ONLY : mpime
    USE laxlib_descriptor, ONLY : la_descriptor, descla_init, laxlib_desc_to_intarray
    USE LAXlib
    implicit none
    !
    TYPE(tester_t) :: test
    !
    integer :: n, m
    complex(DP), allocatable :: h(:,:)            !< full H matrix
    complex(DP), allocatable :: h_save(:,:)       !< full H, used to check consistence across calls
    real(DP), allocatable    :: e(:)              !< full set of eigenvalues
    complex(DP), allocatable :: v(:,:)            !< full eigenvectors
    real(DP), allocatable    :: e_save(:)         !< full set of eigenvalues, used for checks
    complex(DP), allocatable :: v_save(:,:)       !< full set of eigenvectors, used for checks
    TYPE(la_descriptor)      :: desc
    INTEGER :: idesc(LAX_DESC_SIZE)
    !
    integer                  :: i, j, ii, jj, info
    logical                  :: la_proc
    !
    ! Test with a simple Hermitian matrix
    n = 100
    m = n
    !
    ALLOCATE(h(n,n), e(n), v(n,n))
    ALLOCATE(h_save(n,n), e_save(n), v_save(n,n))
    !
    ! Create a random Hermitian matrix on rank 0 and broadcast to all
    IF (mpime == 0) CALL hermitian(n, h)
    CALL mp_bcast(h, 0, world_comm)
    h_save = h
    !
    ! Solve with serial version first to get reference
    CALL diagh(n, m, h_save, e_save, v_save, me_bgrp, root_bgrp, intra_bgrp_comm)
    !
    ! Now test parallel version
    CALL init_parallel_diag(desc, n)
    !
    CALL laxlib_desc_to_intarray( idesc, desc )
    !
    ! Poison v so an unbroadcast rank is detectable.
    v = (1234.5_DP, 6789.0_DP)
    e = 0.d0
    CALL pdiagh( n, h, e, v, idesc )
    !
    CALL test%assert_close( e(1:m), e_save(1:m) )
    !
    ! Run on every rank; |<v_i, v_save_i>| = 1 catches v left
    ! uninitialized on ranks outside the ortho pool.
    DO i = 1, m
       CALL test%assert_close( ABS(SUM(CONJG(v(:,i)) * v_save(:,i))), 1.0_DP )
    END DO
    !
    DEALLOCATE(h, e, v, h_save, e_save, v_save)
    !
  END SUBROUTINE parallel_complex_1
  !
  SUBROUTINE init_parallel_diag(desc, n)

      USE mp_world, ONLY : mpime, nproc, world_comm
      USE laxlib_processors_grid, ONLY : ortho_parent_comm, &
           np_ortho_ => np_ortho, me_ortho_ => me_ortho, &
           ortho_comm_ => ortho_comm, ortho_cntx_ => ortho_cntx
      USE laxlib_descriptor, ONLY : la_descriptor, descla_init
      USE LAXlib
      implicit none
      !
      TYPE(la_descriptor), INTENT(OUT) :: desc
      INTEGER, INTENT(IN) :: n    ! global dimension of the matrix
      !
      INTEGER :: ierr = 0
      INTEGER :: color, key
      !
      INTEGER :: np_ortho(2) = 1  ! size of the processor grid used in ortho
      INTEGER :: me_ortho(2) = 0  ! coordinates of the processors
      INTEGER :: me_ortho1   = 0  ! task id for the ortho group
      INTEGER :: nproc_ortho = 1  ! size of the ortho group:
      INTEGER :: ortho_comm  = 0  ! communicator for the ortho group
      INTEGER :: ortho_cntx  = 0  ! BLACS context for ortho group
      INTEGER :: ortho_row_comm  = 0  ! communicator for the ortho row group
      INTEGER :: ortho_col_comm  = 0  ! communicator for the ortho col group
      INTEGER :: ortho_comm_id = 0 ! id of the ortho_comm
      !
      ortho_parent_comm = world_comm
      !
#if defined __MPI
      !
      CALL grid2d_dims( 'S', nproc, np_ortho(1), np_ortho(2) )
      !
      nproc_ortho = np_ortho(1) * np_ortho(2)
      !
      !  here we choose the first "nproc_ortho" processors
      !
      color = 0
      IF( mpime < nproc_ortho ) color = 1
      !
      key = mpime
      !
      !  initialize the communicator for the new group by splitting the input
      !  communicator
      !
      CALL mpi_comm_split( MPI_COMM_WORLD , color, key, ortho_comm, ierr )
      !
      ! Computes coordinates of the processors, in row maior order
      !
      CALL mpi_comm_rank( ortho_comm, me_ortho1, ierr)
      !
      IF( mpime == 0 .AND. me_ortho1 /= 0 ) &
           CALL lax_error__( " init_ortho_group ", " wrong root task in ortho group ", ierr )
      !
      if( color == 1 ) then
         ! this task belong to the ortho_group compute its coordinates
         ortho_comm_id = 1
         CALL GRID2D_COORDS( 'R', me_ortho1, np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2) )
         CALL GRID2D_RANK( 'R', np_ortho(1), np_ortho(2), me_ortho(1), me_ortho(2), ierr )
         IF( ierr /= me_ortho1 ) &
              CALL lax_error__( " init_ortho_group ", " wrong task coordinates in ortho group ", ierr )
         IF( me_ortho1 /= mpime ) &
              CALL lax_error__( " init_ortho_group ", " wrong rank assignment in ortho group ", ierr )
         CALL mpi_comm_split( ortho_comm , me_ortho(2), me_ortho(1), ortho_col_comm, ierr )
         CALL mpi_comm_split( ortho_comm , me_ortho(1), me_ortho(2), ortho_row_comm, ierr )
         !
         ! Initialize BLACS context for ScaLAPACK
         !
#if defined(__SCALAPACK)
         ortho_cntx = ortho_comm
         call BLACS_GRIDINIT(ortho_cntx, 'R', np_ortho(1), np_ortho(2))
#endif
      else
         ! this task does NOT belong to the ortho_group set dummy values
         ortho_comm_id = 0
         me_ortho(1) = me_ortho1
         me_ortho(2) = me_ortho1
         ortho_cntx = -1
      endif
      !
      ! Set module-level variables for use by laxlib_prdiagh
      !
      np_ortho_ = np_ortho
      me_ortho_ = me_ortho
      ortho_comm_ = ortho_comm
      ortho_cntx_ = ortho_cntx
#else
      ortho_comm_id = 1
#endif
      CALL descla_init( desc, n, n, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )

  END SUBROUTINE init_parallel_diag
  !
end program test_diagh_4
