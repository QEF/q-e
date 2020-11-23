###########################################################
# QE::FOX
###########################################################
if(FOX_ROOT)
    add_library(qe_fox INTERFACE)
    add_library(QE::FOX ALIAS qe_fox)
    qe_install_targets(qe_fox) 
    target_link_libraries(qe_fox INTERFACE "-L${FOX_ROOT}/lib;-lFoX_dom;-lFoX_sax;-lFoX_wxml;-lFoX_common;-lFoX_utils;-lFoX_fsys")
    target_include_directories(qe_fox INTERFACE ${FOX_ROOT}/finclude)
else()
    message(STATUS "Installing QE::FOX via submodule")
    set(fox_targets
        FoX_fsys
        FoX_utils
        FoX_common
        FoX_dom
        FoX_sax
        FoX_wxml)
    set(FoX_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
    qe_git_submodule_update(external/fox)
    add_subdirectory(fox EXCLUDE_FROM_ALL)
    add_library(qe_fox INTERFACE)
    add_library(QE::FOX ALIAS qe_fox)
    target_link_libraries(qe_fox INTERFACE ${fox_targets})
    qe_fix_fortran_modules(${fox_targets})
    qe_install_targets(qe_fox ${fox_targets})
endif()

###########################################################
# QE::libbeef
###########################################################
add_subdirectory(libbeef)
target_link_libraries(qe_libbeef PRIVATE QE::LAPACK)

###########################################################
# QE::WANNIER90
###########################################################
# TODO look for an externally-provided wannier90
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
add_library(QE::WANNIER90 ALIAS qe_wannier90)
target_link_libraries(qe_wannier90
    PRIVATE
        QE::LAPACK)

###########################################################
# wannier_prog.x
###########################################################
set(sources wannier90/src/wannier_prog.F90)
qe_add_executable(qe_wannierprog_exe ${sources})
set_target_properties(qe_wannierprog_exe PROPERTIES OUTPUT_NAME wannier_prog.x)
target_link_libraries(qe_wannierprog_exe
    PRIVATE
        QE::WANNIER90)

###########################################################
# w90chk2chk.x
###########################################################
set(sources wannier90/src/w90chk2chk.F90)
qe_add_executable(qe_w90chk2chk_exe ${sources})
set_target_properties(qe_w90chk2chk_exe PROPERTIES OUTPUT_NAME w90chk2chk.x)
target_link_libraries(qe_w90chk2chk_exe
    PRIVATE
        QE::WANNIER90)

###########################################################
# postw90.x
###########################################################
set(sources wannier90/src/postw90/postw90.F90)
qe_add_executable(qe_wannier90_postw90_exe ${sources})
set_target_properties(qe_wannier90_postw90_exe PROPERTIES OUTPUT_NAME postw90.x)
target_link_libraries(qe_wannier90_postw90_exe
    PRIVATE
        QE::WANNIER90)

###########################################################

qe_install_targets(
    # Libraries
    qe_wannier90
    # Executables
    qe_wannierprog_exe
    qe_w90chk2chk_exe
    qe_wannier90_postw90_exe)

###########################################################
# QE::MBD
###########################################################
if(MBD_ROOT)
    add_library(qe_mbd INTERFACE)
    add_library(QE::MBD ALIAS qe_mbd)
    qe_install_targets(qe_mbd) 
    target_link_libraries(qe_mbd INTERFACE "-L${MBD_ROOT}/lib;-lmbd")
    target_include_directories(qe_mbd INTERFACE ${MBD_ROOT}/finclude)
else()
    message(STATUS "Installing QE::MBD via submodule")
    qe_git_submodule_update(external/mbd)
    if(NOT BUILD_SHARED_LIBS)
        set(BUILD_SHARED_LIBS OFF)
        set(FORCE_BUILD_STATIC_LIBS ON)
    endif()
    add_subdirectory(mbd EXCLUDE_FROM_ALL)
    if(FORCE_BUILD_STATIC_LIBS)
        unset(BUILD_SHARED_LIBS)
    endif()
    add_library(qe_mbd INTERFACE)
    add_library(QE::MBD ALIAS qe_mbd)
    target_link_libraries(qe_mbd INTERFACE Mbd)
    qe_fix_fortran_modules(Mbd)
    qe_install_targets(qe_mbd Mbd)
endif()
