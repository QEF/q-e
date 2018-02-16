PROGRAM test_mp_bcast_iv
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
    INTEGER :: iv(10)
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)

    iv(:) = mpime
    CALL mp_bcast(iv, root, world_comm)
    !
    CALL test%assert_equal(ALL(iv .eq. 0) , .true. )
    !
    iv(:) = mpime
    CALL mp_bcast(iv, nproc-1, world_comm)
    !
    CALL test%assert_equal(ALL(iv .eq. nproc-1) , .true. )
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_iv
