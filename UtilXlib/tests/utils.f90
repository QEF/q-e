SUBROUTINE print_results(test)
#if defined(__MPI)
    USE mpi
#endif
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: itoterr, ierr, me
    !
#if defined(__MPI)
    !
    CALL MPI_REDUCE(test%n_errors, itoterr, 1, MPI_INTEGER, MPI_SUM, &
                        0, MPI_COMM_WORLD, ierr)
    CALL test%assert_equal(0, ierr, fail=.true.)
    CALL test%assert_equal(0, itoterr, fail=.true.)
    
    CALL MPI_Comm_rank(MPI_COMM_WORLD, me, ierr);
    CALL test%assert_equal(0, ierr, fail=.true.)
    
    IF (me .eq. 0) CALL test%print()
#else
    CALL test%print()
#endif
END SUBROUTINE print_results


SUBROUTINE no_test
#if defined(__MPI)
    USE mpi
#endif
    USE tester
    IMPLICIT NONE
    TYPE(tester_t) :: test
    INTEGER :: ierr
    !    
#if defined(__MPI)
    CALL MPI_Init(ierr)
#endif
    CALL test%init()
    CALL print_results(test)
#if defined(__MPI)
    CALL mpi_finalize(ierr)
#endif
END SUBROUTINE no_test
