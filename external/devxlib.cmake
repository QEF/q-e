###########################################################
# DeviceXlib
###########################################################
if(DEVXLIB_ROOT)
    add_library(qe_devxlib INTERFACE)
    # "src" in PATH_SUFFIXES is needed until devxlib supports make install
    find_library(DEVXLIB_LIB
        NAMES devXlib
        HINTS ${DEVXLIB_ROOT}
        PATH_SUFFIXES "lib" "src")
    find_path(DEVXLIB_MOD_PATH
        NAMES device_fbuff_m.mod
        HINTS ${DEVXLIB_ROOT}
        PATH_SUFFIXES "include" "mod" "src")
    target_link_libraries(qe_devxlib INTERFACE ${DEVXLIB_LIB})
    target_include_directories(qe_devxlib INTERFACE ${DEVXLIB_MOD_PATH})

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args("DeviceXlib" REQUIRED_VARS DEVXLIB_LIB DEVXLIB_MOD_PATH)
else()
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
endif()
qe_install_targets(qe_devxlib)
