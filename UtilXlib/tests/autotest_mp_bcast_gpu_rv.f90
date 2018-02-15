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
    ! test variable
    REAL(8), DEVICE :: rv_d(10)
    REAL(8) :: rv_h(10)
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    rv_h(10) = mpime
    rv_d(10) = rv_h(10)
    CALL mp_bcast(rv_d, root, world_comm)
    rv_h(10) = rv_d(10)
    !
    CALL test%assert_equal(ALL(rv_h .eq. 0) , .true. , fail=.true.)
    !
    rv_h(10) = mpime
    rv_d(10) = rv_h(10)
    CALL mp_bcast(rv_d, nproc-1, world_comm)
    rv_h(10) = rv_d(10)
    !
    CALL test%assert_equal(ALL(rv_h .eq. nproc-1) , .true. , fail=.true.)
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
