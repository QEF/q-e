###########################################################
# DeviceXlib
###########################################################
if(QE_DEVICEXLIB_INTERNAL)
    message(STATUS "Installing DeviceXlib via submodule")

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

    target_include_directories(qe_devxlib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/devxlib/src")
    target_include_directories(qe_devxlib PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/devxlib/include")   
else()
    add_library(qe_devxlib INTERFACE)
    find_package(DeviceXlib REQUIRED)
    target_link_libraries(qe_devxlib INTERFACE DeviceXlib::DeviceXlib)
endif()
qe_install_targets(qe_devxlib)
