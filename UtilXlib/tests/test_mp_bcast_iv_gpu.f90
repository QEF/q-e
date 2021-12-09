#if defined(__CUDA)
PROGRAM test_mp_bcast_iv_gpu
!
! Simple program to check the functionalities of test_mp_bcast_i1.
!
    USE cudafor
    USE parallel_include
    USE mp, ONLY : mp_bcast
    USE mp_world, ONLY : mp_world_start, mp_world_end, mpime, &
                          root, nproc, world_comm
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    ! test variable
    INTEGER, DEVICE :: iv_d(10)
    INTEGER :: iv_h(10)
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)
    iv_h(:) = mpime
    iv_d(:) = iv_h(:)
    CALL mp_bcast(iv_d, root, world_comm)
    iv_h(:) = iv_d(:)
    !
    CALL test%assert_equal(ALL(iv_h .eq. 0) , .true. )
    !
    iv_h(:) = mpime
    iv_d(:) = iv_h(:)
    CALL mp_bcast(iv_d, nproc-1, world_comm)
    iv_h(:) = iv_d(:)
    !
    CALL test%assert_equal(ALL(iv_h .eq. nproc-1) , .true. )
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
END PROGRAM test_mp_bcast_iv_gpu
#else
PROGRAM test_mp_bcast_iv_gpu
    CALL no_test()
END PROGRAM test_mp_bcast_iv_gpu
#endif
