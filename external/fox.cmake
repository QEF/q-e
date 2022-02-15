###########################################################
# FoX
###########################################################
add_library(qe_fox INTERFACE)
qe_install_targets(qe_fox)
if(FOX_ROOT)
    find_library(
        FOX_DOM_LIB
        NAMES "FoX_dom"
        HINTS ${FOX_ROOT}
        PATH_SUFFIXES "lib"
        REQUIRED)

    find_library(
        FOX_SAX_LIB
        NAMES "FoX_sax"
        HINTS ${FOX_ROOT}
        PATH_SUFFIXES "lib"
        REQUIRED)

    find_library(
        FOX_WXML_LIB
        NAMES "FoX_wxml"
        HINTS ${FOX_ROOT}
        PATH_SUFFIXES "lib"
        REQUIRED)

    find_library(
        FOX_COMMON_LIB
        NAMES "FoX_common"
        HINTS ${FOX_ROOT}
        PATH_SUFFIXES "lib"
        REQUIRED)

    find_library(
        FOX_UTILS_LIB
        NAMES "FoX_utils"
        HINTS ${FOX_ROOT}
        PATH_SUFFIXES "lib"
        REQUIRED)

    find_library(
        FOX_FSYS_LIB
        NAMES "FoX_fsys"
        HINTS ${FOX_ROOT}
        PATH_SUFFIXES "lib"
        REQUIRED)

    find_path(
        FOX_INC
        NAMES "m_common_io.mod"
        HINTS ${FOX_ROOT}
        PATH_SUFFIXES
            "include"
            "finclude"
            "mod"
            "module"
            "modules"
        REQUIRED)

    target_link_libraries(qe_fox
        INTERFACE
            ${FOX_DOM_LIB}
            ${FOX_SAX_LIB}
            ${FOX_WXML_LIB}
            ${FOX_COMMON_LIB}
            ${FOX_UTILS_LIB}
            ${FOX_FSYS_LIB})
    target_include_directories(qe_fox INTERFACE ${FOX_INC})
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
