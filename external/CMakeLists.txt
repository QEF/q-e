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
    target_link_libraries(qe_lapack INTERFACE ${LAPACK_LIBRARIES})
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
