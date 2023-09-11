qe_add_executable(test_qe_lax_lapack_zdotc test_lapack_zdotc.f90)
set_target_properties(test_qe_lax_lapack_zdotc
    PROPERTIES
        OUTPUT_NAME test_qe_lax_lapack_zdotc.x
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)
target_link_libraries(test_qe_lax_lapack_zdotc
    PRIVATE
        qe_openmp_fortran
        qe_mpi_fortran
        qe_lapack)

add_unit_test(test_qe_lapack_zdotc 1 1 $<TARGET_FILE:test_qe_lax_lapack_zdotc>)
