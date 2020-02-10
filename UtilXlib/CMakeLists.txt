set(sources_util
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
    set_mpi_comm_4_solvers.f90
)

set(sources_util_tests
    tests/tester.f90
    tests/utils.f90
    tests/mp_world.f90
)

# TODO rename all sources_util *.f90 -> *.F90
set_source_files_properties(${sources_util} PROPERTIES COMPILE_FLAGS -cpp)
set_source_files_properties(${sources_util_tests} PROPERTIES COMPILE_FLAGS -cpp)

qe_add_library(qe_utilx ${sources_util})
add_library(QE::UTILX ALIAS qe_utilx)
target_link_libraries(qe_utilx
    PRIVATE
        QE::CLIB
        QE::MPI_Fortran)

qe_install_targets(qe_utilx)

qe_add_library(qe_utilx_tests ${sources_util_tests})
add_library(QE::UTILX::TESTS ALIAS qe_utilx_tests)
target_link_libraries(qe_utilx_tests
    PRIVATE
        QE::UTILX
        QE::MPI_Fortran)

qe_install_targets(qe_utilx_tests)

###########################################################
# tests
###########################################################
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
    test_mp_sum_rv_buffer_gpu
)
# TODO rename all sources *.f90 -> *.F90
FOREACH(NAME ${source_names})
    set(sources tests/${NAME}.f90)
    string(REPLACE "_" "" TGT_NAME ${NAME})
    set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
    qe_add_executable(qe_utilx_tests_${TGT_NAME}_exe ${sources})

    set_target_properties(qe_utilx_tests_${TGT_NAME}_exe PROPERTIES OUTPUT_NAME ${NAME}.x)
    target_link_libraries(qe_utilx_tests_${TGT_NAME}_exe
        PRIVATE
            QE::UTILX
            QE::UTILX::TESTS
            QE::Modules)

    add_test(NAME qe_utilx_tests_${TGT_NAME}_exe COMMAND ${NAME}.x)
ENDFOREACH(NAME)