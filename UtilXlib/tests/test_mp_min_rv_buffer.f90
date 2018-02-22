PROGRAM test_mp_min_rv_buffer
!
! Simple program to check the functionalities of test_mp_min_iv
! with buffer implementation.
!

#if defined(__MPI)
    USE MPI
#endif
    USE mp, ONLY : mp_min
    USE mp_world, ONLY : mp_world_start, mp_world_end, mpime, &
                          root, nproc, world_comm
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: world_group = 0
    ! test variable
    REAL(8) :: rv(200001)
    REAL(8) :: valid(200001)
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)

    rv(:) = mpime + 1
    CALL mp_min(rv, world_comm)
    !
    valid(:) = DBLE(1)
    CALL test%assert_equal( rv, valid )
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
END PROGRAM test_mp_min_rv_buffer
