option(QE_ENABLE_VENDOR_DEPS "enable fallback on vendored deps when none is found via find_package()" ON)

function(qe_git_submodule_update PATH)
    find_package(Git)
    # Old versions of git aren't able to run init+update
    # in one go (via 'git submodule update --init'), we need
    # to call one command for each operation:
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule init -- ${PATH}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update -- ${PATH}
                    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR})
endfunction(qe_git_submodule_update)

###########################################################
# QE::LAPACK
###########################################################
find_package(LAPACK QUIET)
if(LAPACK_FOUND)
    add_library(qe_lapack INTERFACE)
    add_library(QE::LAPACK ALIAS qe_lapack)
    target_link_libraries(qe_lapack INTERFACE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    qe_install_targets(qe_lapack)
else(LAPACK_FOUND)
    if(TARGET QE::LAPACK)
        message(STATUS "Using inherited QE::LAPACK target")
    else(TARGET QE::LAPACK)
        if(QE_ENABLE_VENDOR_DEPS)
            message(STATUS "Installing QE::LAPACK via submodule")
            execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init -- lapack
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
            qe_git_submodule_update(external/lapack)
            add_subdirectory(lapack EXCLUDE_FROM_ALL)
            add_library(qe_lapack INTERFACE)
            add_library(QE::LAPACK ALIAS qe_lapack)
            target_link_libraries(qe_lapack INTERFACE lapack)
            qe_fix_fortran_module_libraries(lapack)
            qe_install_targets(qe_lapack lapack)
        else(QE_ENABLE_VENDOR_DEPS)
            # No dep has been found via find_package,
            # call it again with REQUIRED to make it fail
            # explicitly (hoping in some helpful message)
            find_package(LAPACK REQUIRED)
        endif(QE_ENABLE_VENDOR_DEPS)
    endif(TARGET QE::LAPACK)
endif(LAPACK_FOUND)

###########################################################
# QE::FOX
###########################################################
find_package(FoX QUIET)
if(FoX_FOUND)
    add_library(qe_fox INTERFACE)
    add_library(QE::FOX ALIAS qe_fox)
    target_link_libraries(qe_fox INTERFACE ${FoX_LIBRARIES})
    qe_install_targets(qe_fox)
else(FoX_FOUND)
    if(TARGET QE::FOX)
        message(STATUS "Using inherited QE::FOX target")
    else(TARGET QE::FOX)
        if(QE_ENABLE_VENDOR_DEPS)
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
            qe_fix_fortran_module_libraries(${fox_targets})
            qe_install_targets(qe_fox ${fox_targets})
        else(QE_ENABLE_VENDOR_DEPS)
            # No dep has been found via find_package,
            # call it again with REQUIRED to make it fail
            # explicitly (hoping in some helpful message)
            find_package(FoX REQUIRED)
        endif(QE_ENABLE_VENDOR_DEPS)
    endif(TARGET QE::FOX)
endif(FoX_FOUND)

###########################################################
# QE::IOTK
###########################################################
set(sources
    iotk/src/iotk_dat+COMPLEX3_3.f90
    iotk/src/iotk_attr+LOGICAL4_6.f90
    iotk/src/iotk_attr+REAL1_0.f90
    iotk/src/iotk_attr+INTEGER4_3.f90
    iotk/src/iotk_misc.f90
    iotk/src/iotk_dat+COMPLEX4_6.f90
    iotk/src/iotk_dat+CHARACTER1_6.f90
    iotk/src/iotk_dat+CHARACTER1_3.f90
    iotk/src/iotk_dat+REAL4_3.f90
    iotk/src/iotk_attr+LOGICAL1_6.f90
    iotk/src/iotk_tool.f90
    iotk/src/iotk_attr+LOGICAL2_6.f90
    iotk/src/iotk_attr+LOGICAL3_0.f90
    iotk/src/iotk_attr+COMPLEX1_0.f90
    iotk/src/iotk_dat+INTEGER3_3.f90
    iotk/src/iotk_stream_interf.f90
    iotk/src/iotk_xtox_interf.f90
    iotk/src/iotk_attr+REAL2_0.f90
    iotk/src/iotk_attr+COMPLEX2_6.f90
    iotk/src/iotk_attr+REAL1_6.f90
    iotk/src/iotk_attr+INTEGER3_3.f90
    iotk/src/iotk_attr+LOGICAL4_3.f90
    iotk/src/iotk_fmt.f90
    iotk/src/iotk_dat+INTEGER1_0.f90
    iotk/src/iotk_dat+LOGICAL2_6.f90
    iotk/src/iotk_dat+REAL3_3.f90
    iotk/src/iotk_dat+LOGICAL3_3.f90
    iotk/src/iotk_dat+LOGICAL3_6.f90
    iotk/src/iotk_attr+COMPLEX2_3.f90
    iotk/src/iotk.f90
    iotk/src/iotk_attr+INTEGER3_6.f90
    iotk/src/iotk_attr+LOGICAL1_3.f90
    iotk/src/iotk_str.f90
    iotk/src/iotk_dat+INTEGER4_3.f90
    iotk/src/iotk_dat+REAL1_0.f90
    iotk/src/iotk_attr+INTEGER1_3.f90
    iotk/src/iotk_dat+INTEGER4_0.f90
    iotk/src/iotk_copy.f90
    iotk/src/iotk_attr+COMPLEX1_3.f90
    iotk/src/iotk_dat+COMPLEX2_3.f90
    iotk/src/iotk_attr_interf.f90
    iotk/src/iotk_print_kinds.f90
    iotk/src/iotk_files_interf.f90
    iotk/src/iotk_attr.f90
    iotk/src/iotk_dat_interf.f90
    iotk/src/iotk_attr+INTEGER2_0.f90
    iotk/src/iotk_dat+LOGICAL2_3.f90
    iotk/src/iotk_attr+REAL3_6.f90
    iotk/src/iotk_attr+LOGICAL2_0.f90
    iotk/src/iotk_misc_interf.f90
    iotk/src/iotk_stream.f90
    iotk/src/iotk_attr+COMPLEX3_0.f90
    iotk/src/iotk_dat+COMPLEX3_6.f90
    iotk/src/iotk_attr+INTEGER4_6.f90
    iotk/src/iotk_attr+INTEGER3_0.f90
    iotk/src/iotk_dat+LOGICAL4_0.f90
    iotk/src/iotk_attr+COMPLEX1_6.f90
    iotk/src/iotk_unit_interf.f90
    iotk/src/iotk_dat+LOGICAL1_3.f90
    iotk/src/iotk_attr+LOGICAL3_3.f90
    iotk/src/iotk_dat+LOGICAL2_0.f90
    iotk/src/iotk_attr+LOGICAL1_0.f90
    iotk/src/iotk_attr+REAL4_3.f90
    iotk/src/iotk_dat+REAL2_0.f90
    iotk/src/iotk_attr+COMPLEX3_6.f90
    iotk/src/iotk_dat+INTEGER3_0.f90
    iotk/src/iotk_dat+INTEGER3_6.f90
    iotk/src/iotk_scan.f90
    iotk/src/iotk_str_interf.f90
    iotk/src/iotk_write_interf.f90
    iotk/src/iotk_attr+INTEGER2_3.f90
    iotk/src/iotk_dat+COMPLEX2_0.f90
    iotk/src/iotk_attr+INTEGER4_0.f90
    iotk/src/iotk_attr+REAL3_0.f90
    iotk/src/iotk_dat+REAL1_3.f90
    iotk/src/iotk_dat+LOGICAL4_3.f90
    iotk/src/iotk_dat+REAL1_6.f90
    iotk/src/iotk_dat+INTEGER2_6.f90
    iotk/src/iotk_dat+COMPLEX4_3.f90
    iotk/src/iotk_dat+REAL2_3.f90
    iotk/src/iotk_dat+COMPLEX1_0.f90
    iotk/src/iotk_error_interf.f90
    iotk/src/iotk_attr+REAL4_0.f90
    iotk/src/iotk_dat+LOGICAL3_0.f90
    iotk/src/iotk_dat+REAL3_0.f90
    iotk/src/iotk_dat+REAL2_6.f90
    iotk/src/iotk_error.f90
    iotk/src/iotk_dat+REAL3_6.f90
    iotk/src/iotk_dat+COMPLEX4_0.f90
    iotk/src/iotk_attr+REAL3_3.f90
    iotk/src/iotk_xtox.f90
    iotk/src/iotk_files.f90
    iotk/src/iotk_attr+COMPLEX3_3.f90
    iotk/src/iotk_dat+COMPLEX1_6.f90
    iotk/src/iotk_attr+COMPLEX4_0.f90
    iotk/src/iotk_dat+COMPLEX1_3.f90
    iotk/src/iotk_dat+COMPLEX3_0.f90
    iotk/src/iotk_dat+CHARACTER1_0.f90
    iotk/src/iotk_dat+REAL4_0.f90
    iotk/src/iotk_dat+INTEGER1_6.f90
    iotk/src/iotk_dat+COMPLEX2_6.f90
    iotk/src/iotk_attr+INTEGER1_0.f90
    iotk/src/iotk_write.f90
    iotk/src/iotk_dat+INTEGER2_0.f90
    iotk/src/iotk_attr+LOGICAL3_6.f90
    iotk/src/iotk_dat+LOGICAL1_0.f90
    iotk/src/iotk_attr+REAL1_3.f90
    iotk/src/iotk_attr+REAL2_6.f90
    iotk/src/iotk_attr+LOGICAL2_3.f90
    iotk/src/iotk_attr+REAL2_3.f90
    iotk/src/iotk_dat+LOGICAL1_6.f90
    iotk/src/iotk_dat+INTEGER1_3.f90
    iotk/src/iotk_dat.f90
    iotk/src/iotk_attr+COMPLEX4_3.f90
    iotk/src/iotk_dat+REAL4_6.f90
    iotk/src/iotk_tool_interf.f90
    iotk/src/iotk_module.f90
    iotk/src/iotk_attr+CHARACTER1_0.f90
    iotk/src/iotk_fmt_interf.f90
    iotk/src/iotk_dat+INTEGER2_3.f90
    iotk/src/iotk_attr+COMPLEX4_6.f90
    iotk/src/iotk_attr+REAL4_6.f90
    iotk/src/iotk_unit_list.f90
    iotk/src/iotk_attr+INTEGER2_6.f90
    iotk/src/iotk_scan_interf.f90
    iotk/src/iotk_attr+LOGICAL4_0.f90
    iotk/src/iotk_attr+INTEGER1_6.f90
    iotk/src/iotk_dat+INTEGER4_6.f90
    iotk/src/iotk_base.f90
    iotk/src/iotk_attr+COMPLEX2_0.f90
    iotk/src/iotk_dat+LOGICAL4_6.f90
    iotk/src/iotk_unit.f90)

# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_library(qe_iotk ${sources})
# FIXME: lenght lines of IOTK files are too long for Fortran
target_compile_options(qe_iotk PRIVATE "-ffree-line-length-none")
target_include_directories(qe_iotk
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/iotk/include>
    INTERFACE
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/qe>)
add_library(QE::IOTK ALIAS qe_iotk)

###########################################################

qe_install_targets(qe_iotk)

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
    wannier90/src/postw90/postw90.F90
    wannier90/src/postw90/spin.F90
    wannier90/src/postw90/wan_ham.F90
    wannier90/src/sitesym.F90
    wannier90/src/transport.F90
    wannier90/src/utility.F90
    wannier90/src/w90chk2chk.F90
    wannier90/src/wannierise.F90
    wannier90/src/wannier_lib.F90
    wannier90/src/wannier_prog.F90
    wannier90/src/ws_distance.F90)

# TODO rename all sources *.f90 -> *.F90
set_source_files_properties(${sources} PROPERTIES COMPILE_FLAGS -cpp)
qe_add_library(qe_wannier90 ${sources})
add_library(QE::WANNIER90 ALIAS qe_wannier90)
target_link_libraries(qe_wannier90
    PRIVATE
        QE::LAPACK)
###########################################################

qe_install_targets(qe_wannier90)
