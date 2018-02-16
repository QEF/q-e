#if defined(__CUDA)
PROGRAM test_mp_bcast_rv_gpu
!
! Simple program to check the functionalities of test_mp_bcast_i1.
!
    USE cudafor
#if defined(__MPI)
    USE MPI
#endif
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : mp_world_start, mp_world_end, mpime, &
                          root, nproc, world_comm
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    INTEGER, PARAMETER :: datasize = 10
    !
    ! Stuff for Ranum data
    integer, allocatable :: seed(:)
    integer :: i, n
    REAL :: rnd(datasize)
    !
    ! test variable
    REAL(8), DEVICE :: rv_d(datasize)
    REAL(8) :: rv_h(datasize)
    REAL(8) :: aux_h(datasize)
    
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    rv_h = mpime + 1
    rv_d = rv_h
    CALL mp_bcast(rv_d, root, world_comm)
    rv_h = rv_d
    !
    CALL test%assert_equal(ALL(rv_h .eq. 1) , .true. )
    !
    rv_h = mpime
    rv_d = rv_h
    CALL mp_bcast(rv_d, nproc-1, world_comm)
    rv_h = rv_d
    !
    CALL test%assert_equal(ALL(rv_h .eq. nproc-1) , .true. )
    !
    ! Test against CPU implementation
    CALL random_seed(size = n)
    ALLOCATE(seed(n))
    CALL random_seed(get=seed)
    WRITE (*, *) "Random seed: ", seed
    DEALLOCATE(seed)
    !
    DO i = 0, nproc-1
      CALL RANDOM_NUMBER(rnd)
      rv_h = DBLE ( 10.0 * rnd )
      rv_d = rv_h
      CALL mp_bcast(rv_d, i , world_comm)
      CALL mp_bcast(rv_h, i , world_comm)
      aux_h = rv_d
      CALL test%assert_equal(SUM(rv_h) , SUM(aux_h) )
    END DO
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_rv_gpu
#else
PROGRAM test_mp_bcast_rv_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_rv_gpu
#endif
