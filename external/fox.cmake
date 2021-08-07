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
