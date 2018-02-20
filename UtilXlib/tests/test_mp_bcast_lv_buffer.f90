PROGRAM test_mp_bcast_lv_buffer
!
! Simple program to check the functionalities of test_mp_bcast_lv
! with buffer implementation.
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
    LOGICAL :: lv(200001)
    LOGICAL :: valid
    !
    CALL test%init()
    
#if defined(__MPI)    
    world_group = MPI_COMM_WORLD
#endif
    CALL mp_world_start(world_group)

    lv(:) = MOD(mpime, 2) == 1
    CALL mp_bcast(lv, root, world_comm)
    !
    CALL test%assert_equal(ALL(lv) , .false. )
    !
    lv(:) = MOD(mpime, 2) == 1
    CALL mp_bcast(lv, nproc-1, world_comm)
    !
    valid = MOD(nproc-1, 2) == 1
    CALL test%assert_equal(ALL(lv) , valid )
    !
    CALL print_results(test)
    !
    CALL mp_world_end()
    !
END PROGRAM test_mp_bcast_lv_buffer
