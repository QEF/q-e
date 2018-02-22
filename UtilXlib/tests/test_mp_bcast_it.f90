PROGRAM test_mp_bcast_it
!
! Simple program to check the functionalities of test_mp_bcast_i1.
!

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
    INTEGER :: it(10, 10, 10)
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)

    it(:, :, :) = mpime
    CALL mp_bcast(it, root, world_comm)
    !
    CALL test%assert_equal(ALL(it .eq. 0) , .true. )
    !
    it(:, :, :) = mpime
    CALL mp_bcast(it, nproc-1, world_comm)
    !
    CALL test%assert_equal(ALL(it .eq. nproc-1) , .true. )
    !
    CALL collect_results(test)
    !
    CALL mp_world_end()
    !
    IF (mpime .eq. 0) CALL test%print()
    !
END PROGRAM test_mp_bcast_it
