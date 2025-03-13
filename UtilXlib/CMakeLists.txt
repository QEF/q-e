set(src_util
    clocks_handler.f90
    divide.f90
    data_buffer.f90
    error_handler.f90
    find_free_unit.f90
    fletcher32_mod.f90
    mem_counter.f90
    mp.f90
    mp_base.f90
    mp_bands_util.f90
    nvtx_wrapper.f90
    parallel_include.f90
    print_mem.f90
    util_param.f90
    thread_util.f90
    clib_wrappers.f90
    export_gstart_2_solvers.f90
    set_mpi_comm_4_solvers.f90
    # GPU
    mp_base_gpu.f90)
qe_enable_cuda_fortran("${src_util}")
qe_add_library(qe_utilx ${src_util})

set(src_util_c
    cptimer.c
    copy.c
    c_mkdir.c
    eval_infix.c
    fletcher32.c
    md5.c
    md5_from_file.c
    memstat.c
    memusage.c
    ptrace.c)
add_library(qe_utilx_c ${src_util_c})
target_link_libraries(qe_utilx
    PUBLIC
        qe_utilx_c # cclock is exposed by the mytime module in clocks_handler
        qe_ext_prof_tool # to pass link options, PUBLIC is needed.
        qe_mpi_fortran
    PRIVATE
        qe_openmp_fortran)

qe_install_targets(qe_utilx qe_utilx_c)

# Feature checks
check_function_exists(mallinfo HAVE_MALLINFO)
if(HAVE_MALLINFO)
    target_compile_definitions(qe_utilx_c PRIVATE HAVE_MALLINFO)
endif()

set(src_device_lapack
    device_helper.f90)
qe_enable_cuda_fortran("${src_device_lapack}")
qe_add_library(qe_device_lapack ${src_device_lapack})
target_link_libraries(qe_device_lapack
    PRIVATE
        qe_openacc_fortran
        qe_lapack
        qe_openmp_fortran
        qe_mpi_fortran)

qe_install_targets(qe_device_lapack)

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
        test_offload_macros
        test_mp_count_nodes
        test_mp_bcast_i1
        test_mp_bcast_iv
        test_mp_bcast_im
        test_mp_bcast_it
        test_mp_bcast_iv_buffer
        test_mp_bcast_lv_buffer
        test_mp_bcast_rv_buffer
        test_mp_max_iv_buffer
        test_mp_max_rv_buffer
        test_mp_min_iv_buffer
        test_mp_min_rv_buffer
        test_mp_sum_iv_buffer
        test_mp_sum_rv_buffer)

    if(QE_ENABLE_CUDA)
        set(source_names ${source_names}
            test_mp_bcast_i1_gpu
            test_mp_bcast_iv_gpu
            test_mp_bcast_im_gpu
            test_mp_bcast_it_gpu
            test_mp_bcast_i4d_gpu
            test_mp_bcast_r4d_gpu
            test_mp_bcast_c4d_gpu
            test_mp_bcast_c5d_gpu
            test_mp_bcast_r5d_gpu
            test_mp_bcast_c6d_gpu
            test_mp_bcast_iv_buffer_gpu
            test_mp_bcast_lv_buffer_gpu
            test_mp_bcast_rv_buffer_gpu
            test_mp_max_iv_buffer_gpu
            test_mp_max_rv_buffer_gpu
            test_mp_min_iv_buffer_gpu
            test_mp_min_rv_buffer_gpu
            test_mp_sum_iv_buffer_gpu
            test_mp_sum_rv_buffer_gpu)
    endif()

    # TODO rename all sources *.f90 -> *.F90
    FOREACH(NAME ${source_names})
        set(src_test tests/${NAME}.f90)
        qe_enable_cuda_fortran("${src_test}")
        string(REPLACE "test_" "test_qe_utilx_" TGT_NAME ${NAME})
        qe_add_executable(${TGT_NAME} ${src_test})
        set_target_properties(${TGT_NAME} 
            PROPERTIES
                OUTPUT_NAME ${TGT_NAME}.x 
                RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
        target_link_libraries(${TGT_NAME}
            PRIVATE
                qe_utilx
                qe_utilx_tests
                qe_mpi_fortran)
        add_unit_test(${TGT_NAME} 1 1 $<TARGET_FILE:${TGT_NAME}>)
    ENDFOREACH(NAME)
    target_link_libraries(test_qe_utilx_offload_macros
        PRIVATE
            qe_openacc_fortran
            qe_openmp_fortran)
endif()
