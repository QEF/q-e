#if defined(__CUDA)
PROGRAM test_mp_bcast_it_gpu
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
    INTEGER, DEVICE :: it_d(10, 10, 10)
    INTEGER :: it_h(10, 10, 10)
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    it_h(:, :, :) = mpime
    it_d(:, :, :) = it_h(:, :, :)
    CALL mp_bcast(it_d, root, world_comm)
    it_h(:, :, :) = it_d(:, :, :)
    !
    CALL test%assert_equal(ALL(it_h .eq. 0) , .true. )
    !
    it_h(:, :, :) = mpime
    it_d(:, :, :) = it_h(:, :, :)
    CALL mp_bcast(it_d, nproc-1, world_comm)
    it_h(:, :, :) = it_d(:, :, :)
    !
    CALL test%assert_equal(ALL(it_h .eq. nproc-1) , .true. )
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_it_gpu
#else
PROGRAM test_mp_bcast_it_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_it_gpu
#endif
