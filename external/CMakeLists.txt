###########################################################
# QE::FOX
###########################################################
if(FOX_ROOT)
    add_library(qe_fox INTERFACE)
    add_library(QE::FOX ALIAS qe_fox)
    target_link_libraries(qe_fox INTERFACE "${FOX_ROOT}/lib;-lFoX_fsys;-lFoX_utils;-lFoX_common;-lFoX_wxml;-lFoX_sax;-lFoX_dom")
    target_include_directories(qe_fox INTERFACE ${FOX_ROOT}/include)
elseif(QE_ENABLE_VENDOR_DEPS)
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
else()
    message(FATAL_ERROR "FoX not found! Please specify a FoX root (-DFOX_ROOT=/path/to/fox) or enable the FoX vendoring (-DQE_ENABLE_VENDOR_DEPS=ON)!")
endif()
###########################################################
# QE::libbeef
###########################################################
add_subdirectory(libbeef)
target_link_libraries(qe_libbeef PRIVATE QE::LAPACK)
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
    iotk/src/iotk_attr+INTEGER3_6.f90
    iotk/src/iotk_attr+LOGICAL1_3.f90
    iotk/src/iotk_str.f90
    iotk/src/iotk_dat+INTEGER4_3.f90
    iotk/src/iotk_dat+REAL1_0.f90
    iotk/src/iotk_attr+INTEGER1_3.f90
    iotk/src/iotk_dat+INTEGER4_0.f90
    iotk/src/iotk_attr+COMPLEX1_3.f90
    iotk/src/iotk_dat+COMPLEX2_3.f90
    iotk/src/iotk_attr_interf.f90
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

qe_add_library(qe_iotk ${sources})
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(qe_iotk PRIVATE "-ffree-line-length-none")
endif()
target_include_directories(qe_iotk
    PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/iotk/include>
    INTERFACE
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/qe>)
target_link_libraries(qe_iotk
    PRIVATE
        QE::MPI_Fortran)
add_library(QE::IOTK ALIAS qe_iotk)

###########################################################
# iotk.x
###########################################################
set(sources iotk/src/iotk.f90)
qe_add_executable(qe_iotk_exe ${sources})
set_target_properties(qe_iotk_exe PROPERTIES OUTPUT_NAME iotk.x)
target_link_libraries(qe_iotk_exe
    PRIVATE
        QE::MPI_Fortran
        QE::IOTK)

###########################################################
# iotk_copy.x
###########################################################
set(sources iotk/src/iotk_copy.f90)
qe_add_executable(qe_iotkcopy_exe ${sources})
set_target_properties(qe_iotkcopy_exe PROPERTIES OUTPUT_NAME iotk_copy.x)
target_link_libraries(qe_iotkcopy_exe
    PRIVATE
        QE::MPI_Fortran
        QE::IOTK)

###########################################################
# iotk_print_kinds.x
###########################################################
set(sources iotk/src/iotk_print_kinds.f90)
qe_add_executable(qe_iotkprintkinds_exe ${sources})
set_target_properties(qe_iotkprintkinds_exe PROPERTIES OUTPUT_NAME iotk_print_kinds.x)
target_link_libraries(qe_iotkprintkinds_exe
    PRIVATE
        QE::IOTK)

###########################################################

qe_install_targets(
    # Libraries
    qe_iotk
    # Executables
    qe_iotk_exe
    qe_iotkcopy_exe
    qe_iotkprintkinds_exe
)

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
    qe_wannier90_postw90_exe
)
