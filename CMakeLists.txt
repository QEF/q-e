cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
# Policy CMP0048: The project() command manages VERSION variables
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

project(qe
    VERSION 6.4.1
    DESCRIPTION "ESPRESSO: opEn-Source Package for Research in Electronic Structure, Simulation, and Optimization"
    LANGUAGES Fortran C)

# Add our own custom CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

###########################################################
# QE build framework
# Please use the following functions in place of the
# corresponding CMake builtin
###########################################################

function(qe_add_library LIB)
    add_library(${LIB} ${ARGN})
    # All of the following target modifications make
    # sense on non-interfaces only
    get_target_property(tgt_type ${LIB} TYPE)
    if(NOT ${tgt_type} STREQUAL "INTERFACE_LIBRARY")
        # set module path to tgt_binary_dir/mod
        get_target_property(tgt_binary_dir ${LIB} BINARY_DIR)
        set_target_properties(${LIB}
            PROPERTIES
                Fortran_MODULE_DIRECTORY ${tgt_binary_dir}/mod)
        # make module directory available for clients of LIB 
        target_include_directories(${LIB}
            PUBLIC
                $<BUILD_INTERFACE:${tgt_binary_dir}/mod>
            INTERFACE
                $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/qe>)
    endif()
endfunction(qe_add_library)

function(qe_install_target TGT)
    install(TARGETS ${TGT}
        EXPORT qeTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} # Windows needs RUNTIME also for libraries
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qe)
    # Retrieving non-whitelisted properties leads to an hard
    # error, let's skip the following section for interface
    # targets. See here for details:
    # https://gitlab.kitware.com/cmake/cmake/issues/17640
    get_target_property(tgt_type ${TGT} TYPE)
    if(NOT ${tgt_type} STREQUAL "INTERFACE_LIBRARY")
        # If the target generates Fortran modules, make sure
        # to install them as well to a proper location
        get_target_property(tgt_module_dir ${TGT} Fortran_MODULE_DIRECTORY)
        if(tgt_module_dir)
            install(DIRECTORY ${tgt_module_dir}/
                DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/qe)
        endif()
    endif()
endfunction(qe_install_target)

###########################################################
# Build Type
# Ensure that a specific, default build type is set when 
# none has been explicitly set by the user
###########################################################
set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}"
        CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE
        PROPERTY
            STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

###########################################################
# External dependencies
###########################################################

include(external/external.cmake)

# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

###########################################################
# Build Features
# Each one of the following options boils down to a simple
# compile definition, no additional dependencies are needed
###########################################################
# TODO change all ifdefs throughout code base to match
# cmake options
# TODO symbols beginning with '__' followed by a capital
# character are reserved for standard library use (at
# least in C, not sure about Fortran), change all feature
# macros to avoid weird behaviors
option(QE_ENABLE_TRACE
    "enable execution tracing output" OFF)
option(QE_ENABLE_MPI_INPLACE
    "enable inplace MPI calls (ignored when QE_ENABLE_MPI=OFF)" OFF)
option(QE_ENABLE_MPI_MODULE
    "use MPI via Fortran module instead of mpif.h header inclusion" OFF)
option(QE_ENABLE_BARRIER
    "enable global synchronization between execution units" OFF)

if(QE_ENABLE_TRACE)
    add_definitions(-D__TRACE)
endif()
if(QE_ENABLE_MPI_INPLACE)
    add_definitions(-D__USE_INPLACE_MPI)
endif()
if(QE_ENABLE_MPI_MODULE)
    add_definitions(-D__MPI_MODULE)
endif()
if(QE_ENABLE_BARRIER)
    add_definitions(-D__USE_BARRIER)
endif()

###########################################################
# Build Options
###########################################################

option(QE_ENABLE_TEST
    "enable unit tests" ON)
option(QE_ENABLE_MPI
    "enable distributed execution support via MPI" ON)

###########################################################
# Tests
# Any executable target marked as test runner via
# 'add_test()' will be run by the 'test' make target
###########################################################

if(QE_ENABLE_TEST)
    enable_testing()
endif(QE_ENABLE_TEST)

###########################################################
# MPI
# The following targets will be defined:
qe_add_library(qe_mpi_fortran INTERFACE)
add_library(QE::MPI_Fortran ALIAS qe_mpi_fortran)
qe_add_library(qe_mpi_c INTERFACE)
add_library(QE::MPI_C ALIAS qe_mpi_c)
# Install all internal interface targets to allow
# proper targets export
qe_install_target(qe_mpi_fortran)
qe_install_target(qe_mpi_c)
###########################################################

if(QE_ENABLE_MPI)
    set(mpi_definitions
        # legacy macro used throughout qe
        __MPI
        # skip CXX APIs on openmpi, cause trouble to C APIs
        OMPI_SKIP_MPICXX)
    
    find_package(MPI REQUIRED Fortran C)

    target_link_libraries(qe_mpi_fortran 
        INTERFACE
            MPI::MPI_Fortran)
    set_property(TARGET qe_mpi_fortran
        PROPERTY
            INTERFACE_COMPILE_DEFINITIONS ${mpi_definitions})

    target_link_libraries(qe_mpi_c 
        INTERFACE
            MPI::MPI_C)
    set_property(TARGET qe_mpi_c
        PROPERTY
            INTERFACE_COMPILE_DEFINITIONS ${mpi_definitions})
endif(QE_ENABLE_MPI)

###########################################################
# Components
###########################################################

add_subdirectory(UtilXlib)
add_subdirectory(LAXlib)

###########################################################
# Exports
###########################################################

install(EXPORT qeTargets
        FILE qeTargets.cmake
        NAMESPACE qe::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qe)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    qeConfigVersion.cmake
    VERSION ${PACKAGE_VERSION}
    COMPATIBILITY AnyNewerVersion)

configure_file(
    cmake/qeConfig.cmake.in 
    ${CMAKE_CURRENT_BINARY_DIR}/qeConfig.cmake @ONLY)

install(FILES 
            ${CMAKE_CURRENT_BINARY_DIR}/qeConfigVersion.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/qeConfig.cmake
        DESTINATION
            ${CMAKE_INSTALL_LIBDIR}/cmake/qe)