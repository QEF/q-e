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
        qe_lapack
        qe_scalapack
        qe_openmp_fortran
        qe_mpi_fortran
        qe_utilx
        qe_devxlib)
if(QE_ENABLE_CUDA)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI" OR CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC")
        set(CMAKE_REQUIRED_LINK_OPTIONS "${CUDA_FLAG}lib=cusolver;-fortranlibs")
        check_function_exists(cusolverDnZhegvdx cusolverDnZhegvdx_FOUND)
        if (NOT cusolverDnZhegvdx_FOUND)
            unset(cusolverDnZhegvdx_FOUND CACHE)
            message(FATAL_ERROR "The version of CUDAToolkit chosen by the PGI/NVHPC compiler internally"
                                " doesn't contain cusolverDnZhegvdx. cuSOLVER features used by LAXLib are"
                                " only supported since CUDAToolkit 10.1 release. Use a newer compiler"
                                " or select a newer CUDAToolkit internal to the PGI/NVHPC compiler.")
        endif()
    else()
        if(CUDAToolkit_VERSION VERSION_LESS 10.1)
            message(FATAL_ERROR "cuSOLVER for LAXLib is only supported from CUDA compiler 10.1")
        endif()
    endif()
    target_link_libraries(qe_lax
        PRIVATE
            CUDA::cusolver
            CUDA::cublas)
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
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
    target_link_libraries(qe_lax_test
        PRIVATE
            qe_mpi_fortran
            qe_lax
            qe_lapack)
    add_unit_test(test_qe_lax-r1-t1 1 1 $<TARGET_FILE:qe_lax_test>)
    add_unit_test(test_qe_lax-r1-t3 1 3 $<TARGET_FILE:qe_lax_test>)
    add_unit_test(test_qe_lax-r4-t1 4 1 $<TARGET_FILE:qe_lax_test>)
    add_unit_test(test_qe_lax-r9-t2 9 2 $<TARGET_FILE:qe_lax_test>)

    add_subdirectory(tests)
endif(QE_ENABLE_TEST)
