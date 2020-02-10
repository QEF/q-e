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

###########################################################
# test.x
###########################################################
set(sources test_code.f90)
# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_executable(qe_atomic_test_exe ${sources})
set_target_properties(qe_atomic_test_exe PROPERTIES OUTPUT_NAME test.x)
target_link_libraries(qe_atomic_test_exe
    PRIVATE
        QE::DFTD3)
# FIXME: Does it work?
#add_test(NAME test COMMAND test.x)

qe_install_targets(
    qe_dftd3
)