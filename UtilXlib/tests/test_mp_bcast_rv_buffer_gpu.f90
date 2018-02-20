#if defined(__CUDA)
PROGRAM test_mp_bcast_rv_buffer_gpu
!
! Simple program to check the functionalities of test_mp_bcast_rv
! with buffer implementation.
!

#if defined(__MPI)
    USE MPI
#endif
    USE cudafor
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : mp_world_start, mp_world_end, mpime, &
                          root, nproc, world_comm
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    ! test variable
    REAL(8), DEVICE :: rv_d(200001)
    REAL(8) :: rv_h(200001)
    
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)

    rv_h(:) = mpime+0.5
    !
    rv_d = rv_h
    CALL mp_bcast(rv_d, root, world_comm)
    rv_h = rv_d
    !
    CALL test%assert_equal(ALL(rv_h .eq. 0.5) , .true. )
    !
    rv_h(:) = mpime+0.5
    rv_d = rv_h
    CALL mp_bcast(rv_d, nproc-1, world_comm)
    rv_h = rv_d
    !
    CALL test%assert_equal(ALL(rv_h .eq. nproc-0.5) , .true. )
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_rv_buffer_gpu
#else
PROGRAM test_mp_bcast_rv_buffer_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_rv_buffer_gpu
#endif
