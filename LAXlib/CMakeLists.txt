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
    test.f90
    transto.f90
    zhpev_drv.f90
    cdiaghg.f90
)

# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)

qe_add_library(qe_lax ${sources})
add_library(QE::LAX ALIAS qe_lax)
target_link_libraries(qe_lax
    PUBLIC
        QE::MPI_Fortran
        QE::LAPACK
        QE::UTILX)
qe_install_targets(qe_lax)

###########################################################
# Tests
# TODO move all tests to a proper location
###########################################################
if(QE_ENABLE_TEST)
    # TODO rename all sources *.f90 -> *.F90
    set_source_files_properties(test.f90 PROPERTIES COMPILE_FLAGS -cpp)
    qe_add_executable(qe_lax_test test.f90)
    target_link_libraries(qe_lax_test
        PRIVATE
            QE::MPI_Fortran
            QE::LAX)
    add_test(NAME test-lax COMMAND qe_lax_test)
endif(QE_ENABLE_TEST)
