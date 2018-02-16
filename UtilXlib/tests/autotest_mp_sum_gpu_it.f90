#if defined(__CUDA)
PROGRAM test_mp_sum_it_gpu
!
! Simple program to check the functionalities of test_mp_sum_i1.
!
    USE cudafor
#if defined(__MPI)
    USE MPI
#endif
    USE mp, ONLY : mp_sum
    USE mp_world, ONLY : mp_world_start, mp_world_end, mpime, &
                          root, nproc, world_comm
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0, valid_sum, rnk
    INTEGER, PARAMETER :: datasize = 10
    ! test variable
    INTEGER, DEVICE :: it_d(datasize,datasize,datasize)
    INTEGER :: it_h(datasize,datasize,datasize)
    INTEGER :: aux_h(datasize,datasize,datasize)
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    it_h = mpime
    it_d = it_h
    CALL mp_sum(it_d, world_comm)
    it_h = it_d
    !
    ! The sum of n rank is (zero based)
    !  sum = (n-1)*n*0.5
    !
    ! For a rank N matrix is 2^N * sum
    !
    rnk = SIZE(SHAPE(it_h))
    valid_sum = (10**rnk) * (nproc-1)*nproc/2
    !
    CALL test%assert_equal(INT(SUM(it_h )) , valid_sum )
    !
    ! Validate against CPU implementation
    
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_sum_it_gpu
#else
PROGRAM test_mp_sum_it_gpu
    CALL no_test()
END PROGRAM test_mp_sum_it_gpu
#endif
