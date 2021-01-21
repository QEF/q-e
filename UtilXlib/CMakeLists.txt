set(src_util
    clocks_handler.f90
    cuda_util.f90
    divide.f90
    data_buffer.f90
    error_handler.f90
    find_free_unit.f90
    fletcher32_mod.f90
    mem_counter.f90
    mp.f90
    mp_base.f90
    mp_base_gpu.f90
    mp_bands_util.f90
    parallel_include.f90
    util_param.f90
    thread_util.f90
    export_gstart_2_solvers.f90
    set_mpi_comm_4_solvers.f90)
qe_enable_cuda_fortran("${src_util}")

qe_add_library(qe_utilx ${src_util})
target_link_libraries(qe_utilx
    PRIVATE
        qe_openmp_fortran
        qe_clib
        qe_mpi_fortran)

qe_install_targets(qe_utilx)

###########################################################
# tests
###########################################################
if(QE_ENABLE_TEST)
    set(src_util_tests
        tests/tester.f90
        tests/utils.f90
        tests/mp_world.f90)
    qe_enable_cuda_fortran("${src_util_tests}")

    qe_add_library(qe_utilx_tests ${src_util_tests})
    target_link_libraries(qe_utilx_tests
        PRIVATE
            qe_utilx
            qe_mpi_fortran)

    qe_install_targets(qe_utilx_tests)

    set(source_names
        test_mp_count_nodes
        test_mp_bcast_i1
        test_mp_bcast_i1_gpu
        test_mp_bcast_iv
        test_mp_bcast_iv_gpu
        test_mp_bcast_im
        test_mp_bcast_im_gpu
        test_mp_bcast_it
        test_mp_bcast_it_gpu
        test_mp_bcast_i4d_gpu
        test_mp_bcast_r4d_gpu
        test_mp_bcast_c4d_gpu
        test_mp_bcast_c5d_gpu
        test_mp_bcast_r5d_gpu
        test_mp_bcast_c6d_gpu
        test_mp_bcast_iv_buffer
        test_mp_bcast_iv_buffer_gpu
        test_mp_bcast_lv_buffer
        test_mp_bcast_lv_buffer_gpu
        test_mp_bcast_rv_buffer
        test_mp_bcast_rv_buffer_gpu
        test_mp_max_iv_buffer_gpu
        test_mp_max_iv_buffer
        test_mp_max_rv_buffer
        test_mp_max_rv_buffer_gpu
        test_mp_min_iv_buffer_gpu
        test_mp_min_iv_buffer
        test_mp_min_rv_buffer
        test_mp_min_rv_buffer_gpu
        test_mp_sum_iv_buffer
        test_mp_sum_rv_buffer
        test_mp_sum_iv_buffer_gpu
        test_mp_sum_rv_buffer_gpu)

    # TODO rename all sources *.f90 -> *.F90
    set(src_tests)
    FOREACH(NAME ${source_names})
        set(src_test tests/${NAME}.f90)
        set(src_tests ${src_tests} "${src_test}")
        string(REPLACE "test_" "test_qe_utilx_" TGT_NAME ${NAME})
        qe_add_executable(${TGT_NAME} ${src_test})
        set_target_properties(${TGT_NAME} 
            PROPERTIES
                OUTPUT_NAME ${TGT_NAME}.x 
                RUNTIME_OUTPUT_DIRECTORY ${QE_TESTS_DIR})
        target_link_libraries(${TGT_NAME}
            PRIVATE
                qe_utilx
                qe_utilx_tests
                qe_mpi_fortran
                qe_modules)
        add_test(NAME ${TGT_NAME} COMMAND ${QE_TESTS_DIR}/${TGT_NAME}.x)
    ENDFOREACH(NAME)
endif()
