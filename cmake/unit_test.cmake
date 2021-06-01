# Add a unit test
# TESTNAME the name of test intended for ctest filter.
#   Encoding MPI ranks and OpenMP threads are recommended.
#   For example, XXXX-r2-t3
# PROCS number of MPI ranks
# THREADS number of OpenMP threads
# TEST_BINARY the name of test executable
# additional arguments are passed via ${ARGN}
#
function(add_unit_test TESTNAME PROCS THREADS TEST_BINARY)
    if(NOT QE_ENABLE_OPENMP AND THREADS GREATER 1)
        set(TEST_SKIPPED TRUE)
        message(VERBOSE "Disabling test ${TESTNAME} (building without OpenMP)")
    endif()

    if(NOT TEST_SKIPPED AND NOT QE_ENABLE_MPI AND PROCS GREATER 1)
        set(TEST_SKIPPED TRUE)
        message(VERBOSE "Disabling test ${TESTNAME} (building without MPI)")
    endif()

    if(NOT TEST_SKIPPED)
        if(QE_ENABLE_MPI)
            add_test(NAME ${TESTNAME} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${PROCS} ${MPIEXEC_PREFLAGS} ${TEST_BINARY} ${ARGN})
        else()
            add_test(NAME ${TESTNAME} COMMAND ${TEST_BINARY} ${ARGN})
        endif()

        math(EXPR TOT_PROCS "${PROCS} * ${THREADS}")
        set_tests_properties(${TESTNAME} PROPERTIES
                             PROCESSORS ${TOT_PROCS}
                             ENVIRONMENT OMP_NUM_THREADS=${THREADS}
                             PROCESSOR_AFFINITY TRUE
                             LABELS "unit")
    endif()
endfunction()

