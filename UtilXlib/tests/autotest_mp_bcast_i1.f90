PROGRAM test_mp_bcast_i1
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
    INTEGER, PARAMETER :: datasize = 10
    ! test variable
    INTEGER :: i1
    
    !    
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)

    i1 = mpime
    CALL mp_bcast(i1, root, world_comm)
    !
    CALL test%assert_equal((i1 .eq. 0) , .true. )
    !
    i1 = mpime
    CALL mp_bcast(i1, nproc-1, world_comm)
    !
    CALL test%assert_equal((i1 .eq. nproc-1) , .true. )
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_i1
