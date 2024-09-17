set(sources_xclib
    dft_setting_routines.f90
    dft_setting_params.f90
    qe_constants.f90
    qe_dft_list.f90
    qe_drivers_d_gga.f90
    qe_drivers_d_lda_lsda.f90
    qe_drivers_gga.f90
    qe_drivers_lda_lsda.f90
    qe_drivers_mgga.f90
    qe_funct_corr_gga.f90
    qe_funct_corr_lda_lsda.f90
    qe_funct_exch_gga.f90
    qe_funct_exch_lda_lsda.f90
    qe_funct_mgga.f90
    qe_kind.f90
    xclib_error.f90
    xclib_utils_and_para.f90
    xc_beef_interface.f90
    xc_lib.f90
    xc_wrapper_d_gga.f90
    xc_wrapper_d_lda_lsda.f90
    xc_wrapper_gga.f90
    xc_wrapper_lda_lsda.f90
    xc_wrapper_mgga.f90)
    
qe_add_library(qe_xclib ${sources_xclib})
target_link_libraries(qe_xclib
    PUBLIC
        qe_openacc_fortran
        qe_mpi_fortran
    PRIVATE
        qe_openmp_fortran
        qe_external_libxc
        qe_libbeef)

###########################################################
# libbeef target
# Note: Moving it to a subdirectory XClib/libbeef is prefered.
#       But the legacy Makefile is hard to handle.
###########################################################
set(sources_libbeef
    pbecor.c
    beefun.c)
add_library(qe_libbeef ${sources_libbeef})
target_link_libraries(qe_libbeef 
    PUBLIC
        qe_openacc_c
    PRIVATE 
        qe_lapack)

###########################################################
# xc_infos.x
###########################################################
set(sources_xc_infos
qe_dft_refs.f90
xc_infos.f90)
qe_add_executable(qe_xc_infos ${sources_xc_infos})
set_target_properties(qe_xc_infos
    PROPERTIES
    OUTPUT_NAME qe_xc_infos.x
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
target_link_libraries(qe_xc_infos
    PRIVATE
        qe_external_libxc
        qe_xclib)

###########################################################
qe_install_targets(
    # Libraries
    qe_xclib
    qe_libbeef
    # Executables
    qe_xc_infos
)

###########################################################
# xclib_test.x
###########################################################
if(QE_ENABLE_TEST)
    set(sources_xc_ltest
    qe_dft_refs.f90
    xclib_test.f90)
    qe_add_executable(qe_xclib_test ${sources_xc_ltest})
    set_target_properties(qe_xclib_test 
        PROPERTIES
            OUTPUT_NAME qe_xclib_test.x 
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
    target_link_libraries(qe_xclib_test
        PRIVATE
            qe_openmp_fortran
            qe_mpi_fortran
            qe_upflib
            qe_external_libxc
            qe_libbeef
            qe_xclib)

    find_program(BASH_PROGRAM bash)

    if(NOT BASH_PROGRAM)
        message(WARNING "bash not found! Tests under ${CMAKE_CURRENT_SOURCE_DIR} not added.")
    else()
        set(test_input_file_names_no_extension exe_test gen_test test_all_terms test_all_shorts)

        file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/test_input_files/all_terms.xml
             ${CMAKE_CURRENT_BINARY_DIR}/all_terms.xml
             SYMBOLIC COPY_ON_ERROR)
        file(CREATE_LINK ${CMAKE_CURRENT_SOURCE_DIR}/test_input_files/all_shorts.xml
             ${CMAKE_CURRENT_BINARY_DIR}/all_shorts.xml
             SYMBOLIC COPY_ON_ERROR)
        
        foreach(file_name_no_extension IN LISTS test_input_file_names_no_extension)
            add_unit_test(test_qe_xclib_${file_name_no_extension} 1 1 ${BASH_PROGRAM} -c
                          "$<TARGET_FILE:qe_xclib_test> < ${CMAKE_CURRENT_SOURCE_DIR}/test_input_files/${file_name_no_extension}.in")
        endforeach()

        #additional dependency
        set_tests_properties(test_qe_xclib_exe_test PROPERTIES DEPENDS test_qe_xclib_gen_test)
    endif()
endif(QE_ENABLE_TEST)
