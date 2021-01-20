set(sources_xclib
    beefun.c
    pbecor.c
    dft_mod.f90
    qe_constants.f90
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
    xclib_para_m.f90
    xc_beef_interface.f90
    xc_input_params_mod.f90
    xc_lib.f90
    xc_wrapper_d_gga.f90
    xc_wrapper_d_lda_lsda.f90
    xc_wrapper_gga.f90
    xc_wrapper_lda_lsda.f90
    xc_wrapper_mgga.f90)

qe_add_library(qe_xclib ${sources_xclib})
target_link_libraries(qe_xclib
    PRIVATE
        qe_openmp_fortran
        qe_mpi_fortran
        qe_lapack)

qe_install_targets(qe_xclib)

###########################################################
# xclib_test.x
###########################################################
if(QE_ENABLE_TEST)
    qe_add_executable(qe_xclib_test xclib_test.f90)
    set_target_properties(qe_xclib_test 
        PROPERTIES 
            OUTPUT_NAME qe_xclib_test.x 
            RUNTIME_OUTPUT_DIRECTORY ${QE_TESTS_DIR})
    target_link_libraries(qe_xclib_test
        PRIVATE
            qe_openmp_fortran
            qe_mpi_fortran
            qe_lapack
            qe_xclib)
    add_test(NAME test_qe_xclib COMMAND ${QE_TESTS_DIR}/qe_xclib_test.x)
endif(QE_ENABLE_TEST)
