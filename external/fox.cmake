###########################################################
# FoX
###########################################################
add_library(qe_fox INTERFACE)
qe_install_targets(qe_fox)
if(FOX_ROOT)
    find_library(
        FOX_LIB_COMMON
        NAMES FoX_common
        HINTS ${FOX_ROOT}
        PATH_SUFFIXES "lib")

    if(NOT FOX_LIB_COMMON)
        message(FATAL_ERROR "Failed in locating FoX_common library file at <FOX_ROOT>/lib")
    endif()

    get_filename_component(FOX_LIB_DIR ${FOX_LIB_COMMON} DIRECTORY)

    find_path(
        FOXLIB_MOD_PATH
        NAMES m_common_io.mod
        HINTS ${FOX_ROOT}
        PATH_SUFFIXES "finclude")

    if(NOT FOXLIB_MOD_PATH)
        message(FATAL_ERROR "Failed in locating m_common_io.mod at <FOX_ROOT>/finclude")
    endif()

    target_link_libraries(
        qe_fox INTERFACE "-L${FOX_LIB_DIR};-lFoX_dom;-lFoX_sax;-lFoX_wxml;-lFoX_common;-lFoX_utils;-lFoX_fsys")
    target_include_directories(qe_fox INTERFACE ${FOXLIB_MOD_PATH})
else()
    message(STATUS "Installing FoX via submodule")
    set(fox_targets FoX_fsys FoX_utils FoX_common FoX_dom FoX_sax FoX_wxml)
    set(FoX_ENABLE_EXAMPLES
        OFF
        CACHE BOOL "" FORCE)
    qe_git_submodule_update(external/fox)
    add_subdirectory(fox EXCLUDE_FROM_ALL)
    target_link_libraries(qe_fox INTERFACE ${fox_targets})
    qe_fix_fortran_modules(${fox_targets})
    qe_install_targets(${fox_targets})
endif()
