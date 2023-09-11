###########################################################
# FoX
###########################################################
add_library(qe_fox INTERFACE)
qe_install_targets(qe_fox)
if(QE_ENABLE_FOX)
    if(QE_FOX_INTERNAL)
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
    else()
        find_package(FoX REQUIRED)
        target_link_libraries(qe_fox INTERFACE FoX::FoX)
    endif()
    target_compile_definitions(qe_fox INTERFACE __fox)
else()
    target_link_libraries(qe_fox INTERFACE qe_xml)
endif()
