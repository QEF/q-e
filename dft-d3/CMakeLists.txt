set(sources
    common.f90
    sizes.f90
    pars.f90
    core.f90
    api.f90
    dftd3_qe.f90)

qe_add_library(qe_dftd3 ${sources})
add_library(QE::DFTD3 ALIAS qe_dftd3)
target_link_libraries(qe_dftd3
    PRIVATE
        QE::UTILX 
        QE::Modules)

###########################################################
# test.x
###########################################################
set(sources test_code.f90)
qe_add_executable(qe_atomic_test_exe ${sources})
set_target_properties(qe_atomic_test_exe PROPERTIES OUTPUT_NAME qe_dftd3_test_atomic.x RUNTIME_OUTPUT_DIRECTORY ${QE_TESTS_DIR})
target_link_libraries(qe_atomic_test_exe
    PRIVATE
        QE::Modules
        QE::DFTD3)
# FIXME: Does it work?
#add_test(NAME test_qe_dftd3_atomic COMMAND ${QE_TESTS_DIR}/qe_dftd3_test_atomic.x)

qe_install_targets(qe_dftd3)
