#if defined(__CUDA)
PROGRAM test_mp_bcast_rm_gpu
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
    REAL(8), DEVICE :: rm_d(10,10)
    REAL(8) :: rm_h(10,10)
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    rm_h(10,10) = mpime
    rm_d(10,10) = rm_h(10,10)
    CALL mp_bcast(rm_d, root, world_comm)
    rm_h(10,10) = rm_d(10,10)
    !
    CALL test%assert_equal(ALL(rm_h .eq. 0) , .true. , fail=.true.)
    !
    rm_h(10,10) = mpime
    rm_d(10,10) = rm_h(10,10)
    CALL mp_bcast(rm_d, nproc-1, world_comm)
    rm_h(10,10) = rm_d(10,10)
    !
    CALL test%assert_equal(ALL(rm_h .eq. nproc-1) , .true. , fail=.true.)
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_rm_gpu
#else
PROGRAM test_mp_bcast_rm_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_rm_gpu
#endif
