#if defined(__CUDA)
PROGRAM test_mp_bcast_rt_gpu
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
    REAL :: rnd(datasize,datasize,datasize)
    !
    ! test variable
    REAL(8), DEVICE :: rt_d(datasize,datasize,datasize)
    REAL(8) :: rt_h(datasize,datasize,datasize)
    REAL(8) :: aux_h(datasize,datasize,datasize)
    
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    rt_h = mpime + 1
    rt_d = rt_h
    CALL mp_bcast(rt_d, root, world_comm)
    rt_h = rt_d
    !
    CALL test%assert_equal(ALL(rt_h .eq. 1) , .true. )
    !
    rt_h = mpime
    rt_d = rt_h
    CALL mp_bcast(rt_d, nproc-1, world_comm)
    rt_h = rt_d
    !
    CALL test%assert_equal(ALL(rt_h .eq. nproc-1) , .true. )
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
      rt_h = DBLE ( 10.0 * rnd )
      rt_d = rt_h
      CALL mp_bcast(rt_d, i , world_comm)
      CALL mp_bcast(rt_h, i , world_comm)
      aux_h = rt_d
      CALL test%assert_equal(SUM(rt_h) , SUM(aux_h) )
    END DO
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_rt_gpu
#else
PROGRAM test_mp_bcast_rt_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_rt_gpu
#endif
