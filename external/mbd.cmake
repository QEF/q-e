###########################################################
# MBD
###########################################################
add_library(qe_mbd INTERFACE)
qe_install_targets(qe_mbd)
if(MBD_ROOT)
    find_library(
        MBD_LIB
        NAMES "mbd"
        HINTS ${MBD_ROOT}
        PATH_SUFFIXES "lib"
        REQUIRED)

    find_path(
        MBD_INC
        NAMES "mbd.mod"
        HINTS ${MBD_ROOT}
        PATH_SUFFIXES 
            "mbd"
            "include"
            "finclude"
            "mod"
            "module"
            "modules"
        REQUIRED)

    target_link_libraries(qe_mbd INTERFACE "${MBD_LIB}")
    target_include_directories(qe_mbd INTERFACE "${MBD_INC}")
else()
    message(STATUS "Installing MBD via submodule")
    qe_git_submodule_update(external/mbd)
    if(NOT BUILD_SHARED_LIBS)
        set(BUILD_SHARED_LIBS OFF)
        set(FORCE_BUILD_STATIC_LIBS ON)
    endif()
    set(BUILD_TESTING OFF)
    add_subdirectory(mbd)
    unset(BUILD_TESTING)
    if(FORCE_BUILD_STATIC_LIBS)
        unset(BUILD_SHARED_LIBS)
    endif()
    target_link_libraries(qe_mbd INTERFACE Mbd)
endif()
