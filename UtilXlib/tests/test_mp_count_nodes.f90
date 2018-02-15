PROGRAM test_mp_count_nodes
!
! Simple program to check the functionalities of mp_count_nodes.
! Only num_nodes and keys are tested
!

#if defined(__MPI)
    USE MPI
#endif
    USE mp, ONLY : mp_count_nodes
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    !
    INTEGER :: me, num_nodes, color, key, group, ierr, itoterr
    ! These are validated variables
    INTEGER :: shmcomm, valid_rank, valid_count, valid_n_nodes , is_rank0
    group = 0
    valid_n_nodes = 1
    
    valid_rank = 0
    
    CALL test%init()
    
#if defined(__MPI)    
    group = MPI_COMM_WORLD
    !
    CALL MPI_INIT(ierr)
#endif
    !
    CALL mp_count_nodes(num_nodes, color, key, group)
    !
#if defined(__MPI)
    ! Use MPI 3 to validate results
    CALL MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, &
                                MPI_INFO_NULL, shmcomm, ierr);
    CALL MPI_Comm_rank(shmcomm, valid_rank, ierr);
    CALL MPI_Comm_size(shmcomm, valid_count, ierr);
    IF (valid_rank == 0) THEN
       is_rank0 = 1
    ELSE
       is_rank0 = 0
    END IF
    CALL MPI_ALLREDUCE(is_rank0, valid_n_nodes, 1, MPI_INTEGER, MPI_SUM, &
                          MPI_COMM_WORLD, ierr)
    CALL MPI_COMM_FREE(shmcomm, ierr)
    !
    ! try to split using colors and keys
    CALL MPI_Comm_split(MPI_COMM_WORLD, color, key, shmcomm, ierr);
    ! Check split went fine
    CALL test%assert_equal(ierr, 0, fail=.true.)
    !
    CALL MPI_COMM_FREE(shmcomm, ierr)
#endif

    CALL test%assert_equal(valid_n_nodes, num_nodes, fail=.true.)
    CALL test%assert_equal(valid_rank, key, fail=.true.)
    !
    CALL print_results(test)
    !
#if defined(__MPI)
    CALL MPI_Finalize(ierr)
#endif
END PROGRAM test_mp_count_nodes
