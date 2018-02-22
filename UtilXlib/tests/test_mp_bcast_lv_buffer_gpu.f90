#if defined(__CUDA)
PROGRAM test_mp_bcast_lv_buffer_gpu
!
! Simple program to check the functionalities of test_mp_bcast_lv
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
    LOGICAL, DEVICE :: lv_d(200001)
    LOGICAL :: lv_h(200001)
    LOGICAL :: valid
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)

    lv_h(:) = MOD(mpime, 2) == 1
    !
    lv_d = lv_h
    CALL mp_bcast(lv_d, root, world_comm)
    lv_h = lv_d
    !
    CALL test%assert_equal(ALL(lv_h) , .false. )
    !
    lv_h(:) = MOD(mpime, 2) == 1
    lv_d = lv_h
    CALL mp_bcast(lv_d, nproc-1, world_comm)
    lv_h = lv_d
    !
    valid = MOD(nproc-1, 2) == 1
    CALL test%assert_equal(ALL(lv_h) , valid)
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
END PROGRAM test_mp_bcast_lv_buffer_gpu
#else
PROGRAM test_mp_bcast_lv_buffer_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_lv_buffer_gpu
#endif
