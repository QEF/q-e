set(src_dftd3
    common.f90
    sizes.f90
    pars.f90
    core.f90
    api.f90
    dftd3_qe.f90)

qe_add_library(qe_dftd3 ${src_dftd3})
target_link_libraries(qe_dftd3
    PUBLIC
        qe_openacc_fortran
    PRIVATE
        qe_modules
        qe_upflib)

###########################################################
# test.x
###########################################################
if(QE_ENABLE_TEST)
    set(src_dftd3_test test_code.f90)
    qe_add_executable(qe_atomic_test_exe ${src_dftd3_test})
    set_target_properties(qe_atomic_test_exe 
        PROPERTIES 
            OUTPUT_NAME qe_dftd3_test_atomic.x 
            RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
    target_link_libraries(qe_atomic_test_exe
        PRIVATE
            qe_modules
            qe_dftd3)
    # FIXME: Does it work?
    #add_test(NAME test_qe_dftd3_atomic COMMAND $<TARGET_FILE:qe_atomic_test_exe>)
endif()

qe_install_targets(qe_dftd3)
