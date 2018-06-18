#if defined(__CUDA)
PROGRAM test_mp_max_iv_buffer_gpu
!
! Simple program to check the functionalities of test_mp_max_iv
! with buffer implementation.
!

#if defined(__MPI)
    USE MPI
#endif
    USE cudafor
    USE util_param, ONLY : DP
    USE mp, ONLY : mp_max
    USE mp_world, ONLY : mp_world_start, mp_world_end, mpime, &
                          root, nproc, world_comm
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    ! test variable
    INTEGER, DEVICE :: iv_d(200001)
    INTEGER :: iv_h(200001)
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)

    iv_h(:) = mpime + 1
    iv_d = iv_h
    CALL mp_max(iv_d, world_comm)
    iv_h = iv_d
    !
    CALL test%assert_equal( ALL(iv_h == nproc) , .true. )
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
END PROGRAM test_mp_max_iv_buffer_gpu
#else
PROGRAM test_mp_max_iv_buffer_gpu
    CALL no_test()
END PROGRAM test_mp_max_iv_buffer_gpu
#endif
