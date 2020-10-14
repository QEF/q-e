set(sources
    qe_ext_fun.f90
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
    xc_beef_interface.f90
    xc_get_dft.f90
    xc_input_params_mod.f90
    xc_interfaces.f90
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

qe_install_targets(qe_xclib)


