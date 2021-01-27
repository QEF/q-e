
###########################################################
# FoX
###########################################################
if(FOX_ROOT)
    add_library(qe_fox INTERFACE)
    qe_install_targets(qe_fox) 
    target_link_libraries(qe_fox INTERFACE "-L${FOX_ROOT}/lib;-lFoX_dom;-lFoX_sax;-lFoX_wxml;-lFoX_common;-lFoX_utils;-lFoX_fsys")
    target_include_directories(qe_fox INTERFACE ${FOX_ROOT}/include)
else()
    message(STATUS "Installing FoX via submodule")
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
    target_link_libraries(qe_fox INTERFACE ${fox_targets})
    qe_fix_fortran_modules(${fox_targets})
    qe_install_targets(qe_fox ${fox_targets})
endif()

###########################################################
# WANNIER90
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
target_link_libraries(qe_wannier90
    PRIVATE
        qe_lapack)

###########################################################
# wannier_prog.x
###########################################################
set(sources wannier90/src/wannier_prog.F90)
qe_add_executable(qe_wannierprog_exe ${sources})
set_target_properties(qe_wannierprog_exe PROPERTIES OUTPUT_NAME wannier_prog.x)
target_link_libraries(qe_wannierprog_exe
    PRIVATE
        qe_wannier90)

###########################################################
# w90chk2chk.x
###########################################################
set(sources wannier90/src/w90chk2chk.F90)
qe_add_executable(qe_w90chk2chk_exe ${sources})
set_target_properties(qe_w90chk2chk_exe PROPERTIES OUTPUT_NAME w90chk2chk.x)
target_link_libraries(qe_w90chk2chk_exe
    PRIVATE
        qe_wannier90)

###########################################################
# postw90.x
###########################################################
set(sources wannier90/src/postw90/postw90.F90)
qe_add_executable(qe_wannier90_postw90_exe ${sources})
set_target_properties(qe_wannier90_postw90_exe PROPERTIES OUTPUT_NAME postw90.x)
target_link_libraries(qe_wannier90_postw90_exe
    PRIVATE
        qe_wannier90)

###########################################################

qe_install_targets(
    # Libraries
    qe_wannier90
    # Executables
    qe_wannierprog_exe
    qe_w90chk2chk_exe
    qe_wannier90_postw90_exe)

###########################################################
# MBD
###########################################################
if(MBD_ROOT)
    add_library(qe_mbd INTERFACE)
    qe_install_targets(qe_mbd) 
    target_link_libraries(qe_mbd INTERFACE "-L${MBD_ROOT}/lib;-lmbd")
    target_include_directories(qe_mbd INTERFACE ${MBD_ROOT}/include)
else()
    message(STATUS "Installing MBD via submodule")
    qe_git_submodule_update(external/mbd)
    if(NOT BUILD_SHARED_LIBS)
        set(BUILD_SHARED_LIBS OFF)
        set(FORCE_BUILD_STATIC_LIBS ON)
    endif()
    set(BUILD_TESTING OFF)
    add_subdirectory(mbd EXCLUDE_FROM_ALL)
    unset(BUILD_TESTING)
    if(FORCE_BUILD_STATIC_LIBS)
        unset(BUILD_SHARED_LIBS)
    endif()
    add_library(qe_mbd INTERFACE)
    target_link_libraries(qe_mbd INTERFACE Mbd)
    qe_fix_fortran_modules(Mbd)
    qe_install_targets(qe_mbd Mbd)
endif()

###########################################################
# Eigensolver_GPU
###########################################################
if(QE_ENABLE_CUDA)
    qe_git_submodule_update(external/eigensolver_gpu)

    set(src_eigensolver_gpu
        eigensolver_gpu/lib_eigsolve/dsyevd_gpu.F90
        eigensolver_gpu/lib_eigsolve/dsygst_gpu.F90
        eigensolver_gpu/lib_eigsolve/dsygvdx_gpu.F90
        eigensolver_gpu/lib_eigsolve/dsymv_gpu.F90
        eigensolver_gpu/lib_eigsolve/dsytd2_gpu.F90
        eigensolver_gpu/lib_eigsolve/dsytrd_gpu.F90
        eigensolver_gpu/lib_eigsolve/eigsolve_vars.F90
        eigensolver_gpu/lib_eigsolve/toolbox.F90
        eigensolver_gpu/lib_eigsolve/zheevd_gpu.F90
        eigensolver_gpu/lib_eigsolve/zhegst_gpu.F90
        eigensolver_gpu/lib_eigsolve/zhegvdx_gpu.F90
        eigensolver_gpu/lib_eigsolve/zhemv_gpu.F90
        eigensolver_gpu/lib_eigsolve/zhetd2_gpu.F90
        eigensolver_gpu/lib_eigsolve/zhetrd_gpu.F90)

    # See: https://github.com/NVIDIA/Eigensolver_gpu/blob/master/lib_eigsolve/Makefile
    # Comment: the flags "-O3 -mp -Mlarge_arrays" are inherited
    #          from the global flags of PGI compiler
    set(FLAGS -pgf90libs -Mcuda)
    set(FLAGS2 -pgf90libs -Mcuda,maxregcount:64)
    if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        set(FLAGS3 -pgf90libs -Mcuda=cc35,cc60,nordc,maxregcount:255)
    elseif(CMAKE_CUDA_ARCHITECTURES GREATER 60)
        set(FLAGS3 -pgf90libs -Mcuda=cc60,nordc,maxregcount:255)
    else()
        set(FLAGS3 -pgf90libs -Mcuda=nordc,maxregcount:255)
    endif()
    set_source_files_properties(
        eigensolver_gpu/lib_eigsolve/dsyevd_gpu.F90
        eigensolver_gpu/lib_eigsolve/dsygst_gpu.F90
        eigensolver_gpu/lib_eigsolve/dsygvdx_gpu.F90
        eigensolver_gpu/lib_eigsolve/dsytd2_gpu.F90
        eigensolver_gpu/lib_eigsolve/dsytrd_gpu.F90
        eigensolver_gpu/lib_eigsolve/eigsolve_vars.F90
        eigensolver_gpu/lib_eigsolve/toolbox.F90
        eigensolver_gpu/lib_eigsolve/zheevd_gpu.F90
        eigensolver_gpu/lib_eigsolve/zhegst_gpu.F90
        eigensolver_gpu/lib_eigsolve/zhegvdx_gpu.F90
        eigensolver_gpu/lib_eigsolve/zhetrd_gpu.F90
        PROPERTIES 
            COMPILE_OPTIONS "${FLAGS}")
    set_source_files_properties(
        eigensolver_gpu/lib_eigsolve/zhetd2_gpu.F90
        eigensolver_gpu/lib_eigsolve/dsymv_gpu.F90
        PROPERTIES 
            COMPILE_OPTIONS "${FLAGS2}")
    set_source_files_properties(
        eigensolver_gpu/lib_eigsolve/zhemv_gpu.F90
        PROPERTIES 
            COMPILE_OPTIONS "${FLAGS3}")

    add_library(qe_eigensolver_gpu ${src_eigensolver_gpu})
    qe_fix_fortran_modules(qe_eigensolver_gpu)

    target_link_libraries(qe_eigensolver_gpu
        PRIVATE
            qe_openmp_fortran)
    if(QE_ENABLE_CUDA)
        target_link_libraries(qe_eigensolver_gpu
        PRIVATE 
            CUDA::cusolver)
    endif()

    qe_install_targets(qe_eigensolver_gpu)
endif()

###########################################################
# DeviceXlib
###########################################################
qe_git_submodule_update(external/devxlib)

set(src_devxlib
    devxlib/src/deviceXlib_mod.f90
    devxlib/src/device_memcpy.f90
    devxlib/src/device_memcpy_mod.f90
    devxlib/src/device_auxfunc.f90
    devxlib/src/device_auxfunc_mod.f90
    devxlib/src/device_fbuff.f90
    devxlib/src/device_fbuff_mod.f90
    devxlib/src/timer_mod.f90
    devxlib/src/timer.c)
qe_enable_cuda_fortran("${src_devxlib}")

qe_add_library(qe_devxlib ${src_devxlib})

target_include_directories(qe_devxlib 
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/devxlib/src")
target_include_directories(qe_devxlib 
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/devxlib/include")

qe_install_targets(qe_devxlib)
