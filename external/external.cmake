option(QE_ENABLE_VENDOR_DEPS "enable fallback on vendored deps when none is found via find_package()" ON)

###########################################################
# QE::LAPACK
###########################################################
find_package(LAPACK QUIET)
if(LAPACK_FOUND)
    add_library(qe_lapack INTERFACE)
    add_library(QE::LAPACK ALIAS qe_lapack)
    target_link_libraries(qe_lapack INTERFACE ${LAPACK_LIBRARIES})
else(LAPACK_FOUND)
    if(TARGET QE::LAPACK)
        message(STATUS "Using inherited QE::LAPACK target")
    else(TARGET QE::LAPACK)
        if(QE_ENABLE_VENDOR_DEPS)
            message(STATUS "Installing QE::LAPACK via submodule")
            execute_process(COMMAND git submodule update --init -- external/lapack
                            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
            add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/lapack EXCLUDE_FROM_ALL)
            add_library(QE::LAPACK ALIAS lapack)
        else(QE_ENABLE_VENDOR_DEPS)
            # No dep has been found via find_package,
            # call it again with REQUIRED to make it fail
            # explicitly (hoping in some helpful message)
            find_package(LAPACK REQUIRED)
        endif(QE_ENABLE_VENDOR_DEPS)
    endif(TARGET QE::LAPACK)
endif(LAPACK_FOUND)
