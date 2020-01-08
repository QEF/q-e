set(sources
    common.f90
    sizes.f90
    pars.f90
    core.f90
    api.f90
    dftd3_qe.f90)

# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)

qe_add_library(qe_dftd3 ${sources})
add_library(QE::DFTD3 ALIAS qe_dftd3)
target_link_libraries(qe_dftd3
    PUBLIC
        QE::UTILX 
        QE::Modules)
qe_install_target(qe_dftd3)

