###########################################################
# Tests
###########################################################
set(src_fftx_test fft_test.f90)
qe_add_executable(qe_fftx_test ${src_fftx_test})
set_target_properties(qe_fftx_test 
    PROPERTIES 
        OUTPUT_NAME qe_fftx_test.x 
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
target_link_libraries(qe_fftx_test
    PRIVATE
        qe_openmp_fortran
        qe_mpi_fortran
        qe_fftx)
add_unit_test(test_qe_fftx-r1-t1 1 1 $<TARGET_FILE:qe_fftx_test>)
add_unit_test(test_qe_fftx-r1-t3 1 3 $<TARGET_FILE:qe_fftx_test>)
add_unit_test(test_qe_fftx-r3-t1 3 1 $<TARGET_FILE:qe_fftx_test>)
add_unit_test(test_qe_fftx-r3-t2 3 2 $<TARGET_FILE:qe_fftx_test>)

set(test_common_src tester.f90 utils.f90 sort.f90 recips.f90)
qe_add_library(qe_fftx_test_common ${test_common_src})
target_link_libraries(qe_fftx_test_common PRIVATE qe_fftx qe_mpi_fortran)

set(source_names fft_scalar_gpu fft_scatter_mod_gpu  fwinv_gpu)
foreach(NAME ${source_names})
    set(TEST_SOURCE_FILE test_${NAME}.f90)
    qe_enable_cuda_fortran("${TEST_SOURCE_FILE}")
    qe_add_executable(test_qe_fftx_${NAME} ${TEST_SOURCE_FILE})
    set_target_properties(test_qe_fftx_${NAME}
    PROPERTIES
        OUTPUT_NAME test_qe_fftx_${NAME}.x
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../bin)

    target_link_libraries(test_qe_fftx_${NAME}
        PRIVATE
            qe_openmp_fortran
            qe_mpi_fortran
            qe_fftx_test_common
            qe_fftx)

    add_unit_test(test_qe_fftx_${NAME}-r1-t1 1 1 $<TARGET_FILE:test_qe_fftx_${NAME}>)
    add_unit_test(test_qe_fftx_${NAME}-r1-t3 1 3 $<TARGET_FILE:test_qe_fftx_${NAME}>)
    add_unit_test(test_qe_fftx_${NAME}-r3-t1 3 1 $<TARGET_FILE:test_qe_fftx_${NAME}>)
    add_unit_test(test_qe_fftx_${NAME}-r3-t2 3 2 $<TARGET_FILE:test_qe_fftx_${NAME}>)
endforeach()
