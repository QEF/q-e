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

# Test utils library
set(test_utils_srcs
    tester.f90
    mp_world.f90
    utils.f90
    test_io.f90
    test_helpers.f90)
qe_add_library(qe_lax_test_utils ${test_utils_srcs})
target_link_libraries(qe_lax_test_utils
    PRIVATE
        qe_lax
        qe_utilx
        qe_mpi_fortran)

# Test diagh 1 - basic unit tests
qe_add_executable(test_qe_lax_diagh_1 test_diagh_1.f90)
set_target_properties(test_qe_lax_diagh_1
    PROPERTIES
        OUTPUT_NAME test_diagh_1.x
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)
target_link_libraries(test_qe_lax_diagh_1
    PRIVATE
        qe_lax_test_utils
        qe_lax
        qe_utilx
        qe_mpi_fortran
        qe_lapack)
add_unit_test(test_qe_lax_diagh_1 1 1 $<TARGET_FILE:test_qe_lax_diagh_1>)

# Test diagh 2 - large matrix tests
qe_add_executable(test_qe_lax_diagh_2 test_diagh_2.f90)
set_target_properties(test_qe_lax_diagh_2
    PROPERTIES
        OUTPUT_NAME test_diagh_2.x
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)
target_link_libraries(test_qe_lax_diagh_2
    PRIVATE
        qe_lax_test_utils
        qe_lax
        qe_utilx
        qe_mpi_fortran
        qe_lapack)
add_unit_test(test_qe_lax_diagh_2 1 1 $<TARGET_FILE:test_qe_lax_diagh_2>)

# Test diagh 4 - parallel tests
qe_add_executable(test_qe_lax_diagh_4 test_diagh_4.f90)
set_target_properties(test_qe_lax_diagh_4
    PROPERTIES
        OUTPUT_NAME test_diagh_4.x
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)
target_link_libraries(test_qe_lax_diagh_4
    PRIVATE
        qe_lax_test_utils
        qe_lax
        qe_utilx
        qe_mpi_fortran
        qe_lapack)
add_unit_test(test_qe_lax_diagh_4-r1-t1 1 1 $<TARGET_FILE:test_qe_lax_diagh_4>)
# np=2 exercises the partial-ortho-pool path: GRID2D_DIMS('S', 2, ...) gives
# a 1x1 pool, leaving rank 1 outside. This is the configuration that
# regressed when laxlib_pcdiagh/laxlib_prdiagh failed to broadcast v.
add_unit_test(test_qe_lax_diagh_4-r2-t1 2 1 $<TARGET_FILE:test_qe_lax_diagh_4>)
add_unit_test(test_qe_lax_diagh_4-r4-t1 4 1 $<TARGET_FILE:test_qe_lax_diagh_4>)
