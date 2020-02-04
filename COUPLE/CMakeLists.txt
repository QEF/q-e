set(sources
    src/libpwscf.f90
    src/libcpv.f90
    src/libqemod.f90)

# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)

qe_add_library(qe_couple ${sources})
add_library(QE::COUPLE ALIAS qe_couple)
target_link_libraries(qe_couple
    PUBLIC
        QE::Modules
        QE::LAX
        QE::CPV
        QE::MPI_Fortran)

###########################################################

qe_install_targets(qe_couple)