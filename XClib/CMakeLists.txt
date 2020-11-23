set(sources
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

qe_add_library(qe_xclib ${sources})
add_library(QE::XClib ALIAS qe_xclib)
target_link_libraries(qe_xclib
    PRIVATE
        QE::OpenMP_Fortran
        QE::MPI_Fortran
    PUBLIC
        QE::UTILX)

###########################################################
# xclib_test.x
###########################################################
#set(sources xclib_test.f90)
#qe_add_executable(qe_xclib_xclibtest_exe ${sources})
#set_target_properties(qe_xclib_xclibtest_exe PROPERTIES OUTPUT_NAME xclib_test.x)
#target_link_libraries(qe_xclib_xclibtest_exe
#    PRIVATE
#        QE::XClib)

###########################################################

qe_install_targets(
    # Libraries
    qe_xclib
    # Executables
    # qe_xclib_xclibtest_exe
)

