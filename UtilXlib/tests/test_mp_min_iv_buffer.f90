PROGRAM test_mp_min_iv_buffer
!
! Simple program to check the functionalities of test_mp_min_iv
! with buffer implementation.
!

    USE parallel_include
    USE mp, ONLY : mp_min
    USE mp_world, ONLY : mp_world_start, mp_world_end, mpime, &
                          root, nproc, world_comm
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    ! test variable
    INTEGER :: iv(200001)
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)

    iv(:) = mpime + 1
    CALL mp_min(iv, world_comm)
    !
    CALL test%assert_equal(ALL (iv == 1), .true. )
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
END PROGRAM test_mp_min_iv_buffer
