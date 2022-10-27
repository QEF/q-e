###########################################################
# WANNIER90
###########################################################
if(QE_WANNIER90_INTERNAL)
    message(STATUS "Installing Wannier90 via submodule")
    
    qe_git_submodule_update(external/wannier90)

    set(sources
        wannier90/src/comms.F90
        wannier90/src/constants.F90
        wannier90/src/disentangle.F90
        wannier90/src/hamiltonian.F90
        wannier90/src/io.F90
        wannier90/src/kmesh.F90
        wannier90/src/overlap.F90
        wannier90/src/parameters.F90
        wannier90/src/plot.F90
        wannier90/src/postw90/berry.F90
        wannier90/src/postw90/boltzwann.F90
        wannier90/src/postw90/dos.F90
        wannier90/src/postw90/geninterp.F90
        wannier90/src/postw90/get_oper.F90
        wannier90/src/postw90/gyrotropic.F90
        wannier90/src/postw90/kpath.F90
        wannier90/src/postw90/kslice.F90
        wannier90/src/postw90/postw90_common.F90
        wannier90/src/postw90/spin.F90
        wannier90/src/postw90/wan_ham.F90
        wannier90/src/sitesym.F90
        wannier90/src/transport.F90
        wannier90/src/utility.F90
        wannier90/src/wannierise.F90
        wannier90/src/wannier_lib.F90
        wannier90/src/ws_distance.F90)

    qe_add_library(qe_wannier90 ${sources})
    target_link_libraries(qe_wannier90 PRIVATE qe_lapack)

    ###########################################################
    # wannier90.x
    ###########################################################
    set(sources wannier90/src/wannier_prog.F90)
    qe_add_executable(qe_wannier90_exe ${sources})
    set_target_properties(qe_wannier90_exe PROPERTIES OUTPUT_NAME wannier90.x)
    target_link_libraries(qe_wannier90_exe PRIVATE qe_wannier90)

    ###########################################################
    # w90chk2chk.x
    ###########################################################
    set(sources wannier90/src/w90chk2chk.F90)
    qe_add_executable(qe_w90chk2chk_exe ${sources})
    set_target_properties(qe_w90chk2chk_exe PROPERTIES OUTPUT_NAME w90chk2chk.x)
    target_link_libraries(qe_w90chk2chk_exe PRIVATE qe_wannier90)

    ###########################################################
    # postw90.x
    ###########################################################
    set(sources wannier90/src/postw90/postw90.F90)
    qe_add_executable(qe_wannier90_postw90_exe ${sources})
    set_target_properties(qe_wannier90_postw90_exe PROPERTIES OUTPUT_NAME postw90.x)
    target_link_libraries(qe_wannier90_postw90_exe PRIVATE qe_wannier90)

    ###########################################################

    add_custom_target(w90
        DEPENDS
            qe_wannier90 qe_wannier90_exe qe_w90chk2chk_exe qe_wannier90_postw90_exe
        COMMENT
            "Maximally localised Wannier Functions")

    qe_install_targets(
        # Libraries
        qe_wannier90
        # Executables
        qe_wannier90_exe qe_w90chk2chk_exe qe_wannier90_postw90_exe)
else()
    add_library(qe_wannier90 INTERFACE)
    qe_install_targets(qe_wannier90)
    find_package(Wannier90 REQUIRED)
    target_link_libraries(qe_wannier90 INTERFACE Wannier90::Wannier90)
endif()
