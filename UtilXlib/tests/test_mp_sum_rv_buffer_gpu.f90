#if defined(__CUDA)
PROGRAM test_mp_sum_rv_buffer_gpu
!
! Simple program to check the functionalities of test_mp_sum_rv
! with buffer implementation.
!

#if defined(__MPI)
    USE MPI
#endif
    USE cudafor
    USE util_param, ONLY : DP
    USE mp, ONLY : mp_sum
    USE mp_world, ONLY : mp_world_start, mp_world_end, mpime, &
                          root, nproc, world_comm
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    ! test variable
    INTEGER, DEVICE :: rv_d(200001)
    REAL(DP) :: rv_h(200001)
    REAL(DP) :: valid(200001)
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)

    rv_h(:) = mpime + 1
    rv_d = rv_h
    CALL mp_sum(rv_d, world_comm)
    rv_h = rv_d
    !
    valid(:) = 0.5*nproc*(nproc+1)
    CALL test%assert_equal(rv_h, valid )
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
END PROGRAM test_mp_sum_rv_buffer_gpu
#else
PROGRAM test_mp_sum_rv_buffer_gpu
    CALL no_test()
END PROGRAM test_mp_sum_rv_buffer_gpu
#endif
