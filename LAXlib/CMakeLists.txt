set(src_lax
    distools.f90
    dspev_drv.f90
    la_error.f90
    la_helper.f90
    la_param.f90
    la_types.f90
    mp_diag.f90
    ptoolkit.f90
    rdiaghg.f90
    transto.f90
    zhpev_drv.f90
    cdiaghg.f90)
qe_enable_cuda_fortran("${src_lax}")

qe_add_library(qe_lax ${src_lax})
target_link_libraries(qe_lax
    PRIVATE
        qe_elpa
        qe_scalapack
        qe_openmp_fortran
        qe_mpi_fortran
        qe_lapack
        qe_utilx
        qe_devxlib)
if(QE_ENABLE_CUDA)
    if(QE_ENABLE_LAXLIB_CUSOLVER)
        if(CUDAToolkit_VERSION VERSION_LESS 11)
            message(FATAL_ERROR "cuSOLVER for LAXLib is only supported from CUDA compiler 11, disable it with '-DQE_ENABLE_LAXLIB_CUSOLVER=OFF'")
        endif()
        target_link_libraries(qe_lax
            PRIVATE
                CUDA::cusolver)
        set_property(TARGET qe_lax APPEND
            PROPERTY
                COMPILE_DEFINITIONS __USE_CUSOLVER)
    else()
        target_link_libraries(qe_lax
            PRIVATE
                qe_eigensolver_gpu)
    endif()
endif()

# LAX relies on header files that should be preprocessed
# and renamed *.h -> *.fh due to the fact that those
# headers use preprocessor directives but they are included 
# via fortran 'include' statements that are managed by the
# compiler frontend long after the actual C preprocessor has run.
set(in_headers
    laxlib.h
    laxlib_hi.h
    laxlib_kinds.h
    laxlib_low.h
    laxlib_mid.h
    laxlib_param.h)
set(lax_include_dir ${CMAKE_CURRENT_BINARY_DIR}/include)
foreach(in_h ${in_headers})
    get_filename_component(in_h_basename ${in_h} NAME_WE)
    set(out_h "${lax_include_dir}/${in_h_basename}.fh")
    qe_preprocess_source(${CMAKE_CURRENT_SOURCE_DIR}/${in_h} ${out_h})
    list(APPEND out_headers ${out_h})
endforeach()
add_custom_target(qe_lax_headers
    DEPENDS ${out_headers}
    VERBATIM)
target_include_directories(qe_lax
    PUBLIC
        $<BUILD_INTERFACE:${lax_include_dir}>
        $<INSTALL_INTERFACE:include/qe>)
set_target_properties(qe_lax
    PROPERTIES PUBLIC_HEADER "${out_headers}")
add_dependencies(qe_lax qe_lax_headers)

qe_install_targets(qe_lax)

###########################################################
# Tests
# TODO move all tests to a proper location
###########################################################
if(QE_ENABLE_TEST)
    set(src_lax_test test.f90)
    qe_enable_cuda_fortran("${src_lax_test}")
    qe_add_executable(qe_lax_test ${src_lax_test})
    set_target_properties(qe_lax_test 
        PROPERTIES 
            OUTPUT_NAME qe_lax_test.x 
            RUNTIME_OUTPUT_DIRECTORY ${QE_TESTS_DIR})
    target_link_libraries(qe_lax_test
        PRIVATE
            qe_mpi_fortran
            qe_lapack
            qe_lax)
    add_test(NAME test_qe_lax 
        COMMAND ${QE_TESTS_DIR}/qe_lax_test.x)
endif(QE_ENABLE_TEST)
