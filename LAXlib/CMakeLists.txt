set(sources
    la_types.f90
    la_error.f90
    cdiaghg.f90
    rdiaghg.f90
    dspev_drv.f90
    ptoolkit.f90
    transto.f90
    distools.f90
    zhpev_drv.f90
    mp_diag.f90
    set_mpi_comm_4_solvers.f90
    export_gstart_2_solvers.f90
    la_param.f90
)

# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)

add_fortran_library(qe_laxlib ${sources})
add_library(QE::LAXLib ALIAS qe_laxlib)
target_link_libraries(qe_laxlib
    PUBLIC
        QE::Features
        QE::MPI_Fortran
        QE::LAPACK
        QE::UtilXLib)

###########################################################
# Tests
# TODO move all tests to a proper location
###########################################################
if(QE_ENABLE_TEST)
    # TODO rename all sources *.f90 -> *.F90
    set_source_files_properties(test.f90 PROPERTIES COMPILE_FLAGS -cpp)
    add_executable(qe_laxlib_test test.f90)
    target_link_libraries(qe_laxlib_test
        PRIVATE
            QE::Features
            QE::MPI_Fortran
            QE::LAXLib)
    add_test(NAME test-laxlib COMMAND qe_laxlib_test)
endif(QE_ENABLE_TEST)