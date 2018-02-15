#if defined(__CUDA)
PROGRAM test_mp_bcast_i1_gpu
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
    ! test variable
    INTEGER, DEVICE :: i_d
    INTEGER :: i_h
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    i_h = mpime
    i_d = i_h
    CALL mp_bcast(i_d, root, world_comm)
    i_h = i_d
    !
    CALL test%assert_equal(i_h, 0 , fail=.true.)
    !
    i_h = mpime
    i_d = i_h
    CALL mp_bcast(i_d, nproc-1, world_comm)
    i_h = i_d
    !
    CALL test%assert_equal(i_h, nproc-1, fail=.true.)
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_i1_gpu
#else
PROGRAM test_mp_bcast_i1_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_i1_gpu
#endif
