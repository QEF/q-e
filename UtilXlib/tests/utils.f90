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

SUBROUTINE save_random_seed(test_name, mpime)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: test_name
    INTEGER, INTENT(IN) :: mpime
    !
    INTEGER, PARAMETER :: in_unit=20, out_unit=21
    CHARACTER(len=80) :: fname
    INTEGER :: n, istat
    INTEGER, ALLOCATABLE :: seed(:)
    !
    CALL random_seed(size = n)
    ALLOCATE(seed(n))

   ! First try if the OS provides a random number generator
    OPEN(UNIT=in_unit, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
    
    IF (istat == 0) THEN
        READ(in_unit) seed
        CLOSE(in_unit)
    ELSE
        ! Fallback to stupid algorithm. Actually we do not really need
        !  high-quality random numbers
        CALL random_seed(get=seed)
        seed = seed + mpime
    END IF
    !
    CALL random_seed(put=seed)
    !
    WRITE(fname, '("rnd_seed_",A,I4.4)') TRIM(test_name), mpime
    fname = TRIM(fname)
    !
    OPEN (UNIT=out_unit,FILE=fname,ACTION="write",STATUS="replace")
    !
    WRITE (out_unit,*) n
    WRITE (out_unit,*) seed
    CLOSE (out_unit)
    DEALLOCATE(seed)
    !
END SUBROUTINE save_random_seed


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
