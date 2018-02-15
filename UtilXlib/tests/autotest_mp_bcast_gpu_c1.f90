#if defined(__CUDA)
PROGRAM test_mp_bcast_c1_gpu
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
    COMPLEX(8), DEVICE :: c1_d
    COMPLEX(8) :: c1_h
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    c1_h = mpime
    c1_d = c1_h
    CALL mp_bcast(c1_d, root, world_comm)
    c1_h = c1_d
    !
    CALL test%assert_equal((c1_h .eq. 0) , .true. , fail=.true.)
    !
    c1_h = mpime
    c1_d = c1_h
    CALL mp_bcast(c1_d, nproc-1, world_comm)
    c1_h = c1_d
    !
    CALL test%assert_equal((c1_h .eq. nproc-1) , .true. , fail=.true.)
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_c1_gpu
#else
PROGRAM test_mp_bcast_c1_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_c1_gpu
#endif
