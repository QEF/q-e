set(sources
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

qe_add_library(qe_lax ${sources})
add_library(QE::LAX ALIAS qe_lax)
target_link_libraries(qe_lax
    PRIVATE
        QE::ELPA
        QE::SCALAPACK
        QE::OpenMP_Fortran
        QE::MPI_Fortran
        QE::LAPACK
        QE::UTILX)

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
    qe_add_executable(qe_lax_test test.f90)
    set_target_properties(qe_lax_test PROPERTIES OUTPUT_NAME qe_lax_test.x RUNTIME_OUTPUT_DIRECTORY ${QE_TESTS_DIR})
    target_link_libraries(qe_lax_test
        PRIVATE
            QE::MPI_Fortran
            QE::LAPACK
            QE::LAX)
    add_test(NAME test_qe_lax COMMAND ${QE_TESTS_DIR}/qe_lax_test.x)
endif(QE_ENABLE_TEST)
