! This test uses the internal parallel diagonalization algorithm of LAXlib
! to solve the problems stored in binary files:
!
!  - ZnOG1.bin
!  - ZnOG2.bin
!  - ZnOK1.bin
!  - ZnOK2.bin
!  - SiGeK1.bin
!  - SiGeK2.bin
!
! If the scalacpak or ELPA driver is used, the test is skipped.
!
#if ! defined(__SCALAPACK)
program test_diaghg_4
#if defined(__MPI)
    USE MPI
#endif
    USE mp,            ONLY : mp_bcast
    USE mp_world,      ONLY : mp_world_start, mp_world_end, mpime, &
                              root, nproc, world_comm
    USE mp_bands_util, ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
    USE tester
    IMPLICIT NONE
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
    USE descriptors, ONLY : la_descriptor, descla_init, descla_local_dims
    USE la_param,    ONLY : DP
    USE test_io
    implicit none
    !
    TYPE(tester_t) :: test
    !
    TYPE(la_descriptor) :: desc
    integer :: ldh, n, m
    real(DP), allocatable :: h(:,:), hdst(:,:) !< full and distributed Hpsi
    real(DP), allocatable :: h_save(:,:)       !< full Hpsi, used to check consistence across calls
    real(DP), allocatable :: s(:,:), sdst(:,:) !< full and distributed Spsi
    real(DP), allocatable :: s_save(:,:)       !< full Spsi, used to check consistence across calls
    real(DP), allocatable    :: e(:)           !< full set of eigenvalues
    real(DP), allocatable :: v(:,:), vdst(:,:) !< full and distributed eigenvectors
    real(DP), allocatable    :: e_save(:)      !< full set of eigenvalues, used for checks
    real(DP), allocatable :: v_save(:,:)       !< full set of eigenvectors, used for checks
    !
    character(len=20)        :: inputs(2)
    integer                  :: l, i, j, ii, jj, info, nrdst
    logical                  :: la_proc
    !
    inputs = ["ZnOG1.bin", "ZnOG2.bin"]
    !
    DO l=1, SIZE(inputs)
        !
        CALL read_problem(inputs(l), ldh, n, m, h, s, e, v, info)
        !
        IF (info /= 0) THEN
            IF (mpime == 0) print *, "Test with ", inputs(l), " skipped. Input not found."
            CYCLE
        END IF
        !
        ALLOCATE(h_save, SOURCE=h)
        ALLOCATE(s_save, SOURCE=s)
        ALLOCATE(e_save, SOURCE=e)
        ALLOCATE(v_save, SOURCE=v)
        !
        h_save = h
        s_save = s
        !
        CALL init_parallel_diag(desc, n)
        !
        IF( desc%active_node > 0 ) la_proc = .TRUE.
        nrdst = desc%nrcx
        IF (.not. la_proc) nrdst = 1
        !
        v = (0.d0, 0.d0)
        e = 0.d0
        print *, nrdst, n, m
        ALLOCATE( hdst( nrdst , nrdst ), STAT=info )
        ALLOCATE( sdst( nrdst , nrdst ), STAT=info )
        ALLOCATE( vdst( nrdst , nrdst ), STAT=info )
        !
        IF (la_proc) THEN
            DO j = 1, desc%nc ! number of column in the local block of lambda
                DO i = 1, desc%nr ! number of row in the local block of lambda
                   ii = i + desc%ir - 1 ! globla index of the first row in the local block of lambda
                   jj = j + desc%ic - 1 ! global index of the first column in the local block of lambda
                   hdst(i, j) = h(ii, jj)
                   sdst(i, j) = s(ii, jj)
                END DO
            END DO
        END IF
        !
        CALL pdiaghg( n, hdst, sdst, nrdst, e, vdst, desc, .false. )
        !
        DO j = 1, m
            !CALL test%assert_close( v(1:n, j), v_save(1:n, j))
        END DO
        CALL test%assert_close( e(1:m), e_save(1:m) )
        !
        !
        v = (0.d0, 0.d0)
        e = 0.d0
        CALL pdiaghg( n, hdst, sdst, nrdst, e, vdst, desc, .true. )
        !
        DO j = 1, m
            !CALL test%assert_close( v(1:n, j), v_save(1:n, j))
        END DO
        CALL test%assert_close( e(1:m), e_save(1:m))
        !
        DEALLOCATE(h,s,e,v,h_save,s_save,e_save,v_save, hdst, sdst, vdst)
    END DO
    !
  END SUBROUTINE parallel_real_1
  !
  SUBROUTINE parallel_complex_1(test)
    USE mp_world, ONLY : mpime
    USE descriptors, ONLY : la_descriptor, descla_init, descla_local_dims
    USE LAXlib
    USE la_param, ONLY : DP
    USE test_io
    implicit none
    !
    TYPE(tester_t) :: test
    ! 
    integer :: ldh, n, m
    complex(DP), allocatable :: h(:,:), hdst(:,:) !< full and distributed Hpsi
    complex(DP), allocatable :: h_save(:,:)       !< full Hpsi, used to check consistence across calls
    complex(DP), allocatable :: s(:,:), sdst(:,:) !< full and distributed Spsi
    complex(DP), allocatable :: s_save(:,:)       !< full Spsi, used to check consistence across calls
    real(DP), allocatable    :: e(:)              !< full set of eigenvalues
    complex(DP), allocatable :: v(:,:), vdst(:,:) !< full and distributed eigenvectors
    real(DP), allocatable    :: e_save(:)         !< full set of eigenvalues, used for checks
    complex(DP), allocatable :: v_save(:,:)       !< full set of eigenvectors, used for checks
    TYPE(la_descriptor)      :: desc
    !
    character(len=20)        :: inputs(4)
    integer                  :: l, i, j, ii, jj, info, nrdst
    logical                  :: la_proc
    !
    inputs = ["ZnOK1.bin ", &
              "ZnOK2.bin ", &
              "SiGeK1.bin", &
              "SiGeK2.bin"]
    !
    DO l=1, SIZE(inputs)
        !
        CALL read_problem(inputs(l), ldh, n, m, h, s, e, v, info)
        !
        IF (info /= 0) THEN
            IF (mpime == 0) print *, "Test with ", inputs(l), " skipped. Input not found."
            CYCLE
        END IF
        !
        ALLOCATE(h_save, SOURCE=h)
        ALLOCATE(s_save, SOURCE=s)
        ALLOCATE(e_save, SOURCE=e)
        ALLOCATE(v_save, SOURCE=v)
        !
        h_save = h
        s_save = s
        !
        CALL init_parallel_diag(desc, n)
        !
        IF( desc%active_node > 0 ) la_proc = .TRUE.
        nrdst = desc%nrcx
        IF (.not. la_proc) nrdst = 1
        !
        v = (0.d0, 0.d0)
        e = 0.d0
        print *, nrdst, n, m
        ALLOCATE( hdst( nrdst , nrdst ), STAT=info )
        ALLOCATE( sdst( nrdst , nrdst ), STAT=info )
        ALLOCATE( vdst( nrdst , nrdst ), STAT=info )
        !
        IF (la_proc) THEN
            DO j = 1, desc%nc ! number of column in the local block of lambda
                DO i = 1, desc%nr ! number of row in the local block of lambda
                   ii = i + desc%ir - 1 ! globla index of the first row in the local block of lambda
                   jj = j + desc%ic - 1 ! global index of the first column in the local block of lambda
                   hdst(i, j) = h(ii, jj)
                   sdst(i, j) = s(ii, jj)
                END DO
            END DO
        END IF
        !
        CALL pdiaghg( n, hdst, sdst, nrdst, e, vdst, desc, .false. )
        !
        DO j = 1, m
            !CALL test%assert_close( v(1:n, j), v_save(1:n, j))
        END DO
        CALL test%assert_close( e(1:m), e_save(1:m) )
        !
        !
        v = (0.d0, 0.d0)
        e = 0.d0
        CALL pdiaghg( n, hdst, sdst, nrdst, e, vdst, desc, .true. )
        !
        DO j = 1, m
            !CALL test%assert_close( v(1:n, j), v_save(1:n, j))
        END DO
        CALL test%assert_close( e(1:m), e_save(1:m))
        !
        DEALLOCATE(h,s,e,v,h_save,s_save,e_save,v_save, hdst, sdst, vdst)
    END DO
    !
  END SUBROUTINE parallel_complex_1
  !
  SUBROUTINE init_parallel_diag(desc, n)
  
      USE mp_world, ONLY : mpime, nproc, world_comm
      USE mp_diag,  ONLY : ortho_parent_comm
      USE descriptors, ONLY : la_descriptor, descla_init, descla_local_dims
      USE LAXlib
      USE la_param, ONLY : DP
      implicit none
      !
      TYPE(la_descriptor) :: desc
      INTEGER             :: n    ! global dimension of the matrix
      !
      INTEGER :: ierr = 0
      INTEGER :: color, key
      !
      INTEGER :: np_ortho(2) = 1  ! size of the processor grid used in ortho
      INTEGER :: me_ortho(2) = 0  ! coordinates of the processors
      INTEGER :: me_ortho1   = 0  ! task id for the ortho group
      INTEGER :: nproc_ortho = 1  ! size of the ortho group:
      INTEGER :: ortho_comm  = 0  ! communicator for the ortho group
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
      else
         ! this task does NOT belong to the ortho_group set dummy values
         ortho_comm_id = 0
         me_ortho(1) = me_ortho1
         me_ortho(2) = me_ortho1
      endif
#else
      ortho_comm_id = 1
#endif
      CALL descla_init( desc, n, n, np_ortho, me_ortho, ortho_comm, -1, ortho_comm_id )
      
  END SUBROUTINE init_parallel_diag
  
end program test_diaghg_4
#else
program test_diaghg_4
end program test_diaghg_4
#endif
