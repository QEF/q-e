###########################################################
# Environ
###########################################################
if(QE_ENABLE_ENVIRON)
    if(QE_ENABLE_ENVIRON MATCHES "INTERNAL|internal")
        message(FATAL_ERROR "Not yet implemented. Use QE_ENABLE_ENVIRON=EXTERNAL")
    else()
        message(STATUS "Linking external Environ library")
        add_library(qe_environ INTERFACE)
        find_package(Environ REQUIRED)
        target_link_libraries(qe_environ INTERFACE Environ::Environ)
    endif()
    qe_install_targets(qe_environ)
endif()
