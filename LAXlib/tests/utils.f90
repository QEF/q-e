SUBROUTINE collect_results(test)
#if defined(__MPI)
    USE mpi
#endif
    USE tester
    IMPLICIT NONE
    !
    TYPE(tester_t) :: test
    INTEGER :: itottests, itoterr, ierr, me
    !
#if defined(__MPI)
    !
    CALL MPI_REDUCE(test%n_errors, itoterr, 1, MPI_INTEGER, MPI_SUM, &
                        0, MPI_COMM_WORLD, ierr)
    ! Fail in case MPI fails...
    IF (ierr /= 0) CALL test%assert_equal(0, ierr)
    !
    CALL MPI_REDUCE(test%n_tests, itottests, 1, MPI_INTEGER, MPI_SUM, &
                        0, MPI_COMM_WORLD, ierr)
    ! Fail in case MPI fails...
    IF (ierr /= 0) CALL test%assert_equal(0, ierr)
    !
    test%n_tests  = itottests
    test%n_errors = itoterr
    !
    IF (ierr /= 0) CALL test%assert_equal(0, ierr)
    !
    CALL MPI_Comm_rank(MPI_COMM_WORLD, me, ierr);
    !
    IF (ierr /= 0) CALL test%assert_equal(0, ierr)
    !
#endif
END SUBROUTINE collect_results


SUBROUTINE no_test
#if defined(__MPI)
    USE mpi
#endif
    USE tester
    IMPLICIT NONE
    !TYPE(tester_t) :: test
    INTEGER :: ierr
    !    
#if defined(__MPI)
    CALL MPI_Init(ierr)
#endif
    !CALL test%init()
    !CALL print_results(test)
#if defined(__MPI)
    CALL mpi_finalize(ierr)
#endif
END SUBROUTINE no_test
