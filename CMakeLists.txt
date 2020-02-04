cmake_minimum_required(VERSION 3.14 FATAL_ERROR)
# CMake < 3.13 cannot install external targets:
# https://gitlab.kitware.com/cmake/cmake/merge_requests/2152

# CMake < 3.14 generates incorrect dependency graphs with
# alias targets:
# https://gitlab.kitware.com/cmake/cmake/merge_requests/2521

# Policy CMP0048: The project() command manages VERSION variables
set(CMAKE_POLICY_DEFAULT_CMP0048 NEW)

project(qe
    VERSION 6.4.1
    DESCRIPTION "ESPRESSO: opEn-Source Package for Research in Electronic Structure, Simulation, and Optimization"
    LANGUAGES Fortran C)

###########################################################
# Build helpers
###########################################################
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(QE_MANDATORY_TARGETS QE::Definitions CACHE STRING "" FORCE)
include(cmake/qeHelpers.cmake)

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
# Modules
###########################################################
include(CheckFunctionExists)
# Must use GNUInstallDirs to install libraries into correct
# locations on all platforms.
include(GNUInstallDirs)

###########################################################
# Build Options
###########################################################
option(QE_ENABLE_TEST
    "enable unit tests" ON)
option(QE_ENABLE_MPI
    "enable distributed execution support via MPI" ON)

###########################################################
# Build Features
# Each one of the following options boils down to a simple
# compile definition, no additional dependencies are needed
add_library(qe_definitions INTERFACE)
add_library(QE::Definitions ALIAS qe_definitions)
qe_install_targets(qe_definitions)
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

# Disable all configuration headers used to be generated
# by configure (see <qe>/include/)
set_property(TARGET qe_definitions APPEND
    PROPERTY INTERFACE_COMPILE_DEFINITIONS QE_NO_CONFIGURE_HEADERS)

if(QE_ENABLE_TRACE)
    set_property(TARGET qe_definitions APPEND
        PROPERTY INTERFACE_COMPILE_DEFINITIONS __TRACE)
endif()
if(QE_ENABLE_MPI_INPLACE)
    set_property(TARGET qe_definitions APPEND
        PROPERTY INTERFACE_COMPILE_DEFINITIONS __USE_INPLACE_MPI)
endif()
if(QE_ENABLE_MPI_MODULE)
    set_property(TARGET qe_definitions APPEND
        PROPERTY INTERFACE_COMPILE_DEFINITIONS __MPI_MODULE)
endif()
if(QE_ENABLE_BARRIER)
    set_property(TARGET qe_definitions APPEND
        PROPERTY INTERFACE_COMPILE_DEFINITIONS __USE_BARRIER)
endif()

# Feature checks
check_function_exists(mallinfo HAVE_MALLINFO)
if(HAVE_MALLINFO)
    set_property(TARGET qe_definitions APPEND
        PROPERTY INTERFACE_COMPILE_DEFINITIONS HAVE_MALLINFO)
endif()

###########################################################
# MPI
# The following targets will be defined:
add_library(qe_mpi_fortran INTERFACE)
add_library(QE::MPI_Fortran ALIAS qe_mpi_fortran)
add_library(qe_mpi_c INTERFACE)
add_library(QE::MPI_C ALIAS qe_mpi_c)
qe_install_targets(qe_mpi_fortran qe_mpi_c)
###########################################################
if(QE_ENABLE_MPI)
    set(mpi_definitions
        # legacy macro used throughout qe
        __MPI
        # skip CXX APIs on openmpi, cause trouble to C APIs
        OMPI_SKIP_MPICXX)
    
    find_package(MPI REQUIRED Fortran C)

    target_link_libraries(qe_mpi_fortran 
        INTERFACE MPI::MPI_Fortran)
    set_property(TARGET qe_mpi_fortran
        PROPERTY INTERFACE_COMPILE_DEFINITIONS ${mpi_definitions})

    target_link_libraries(qe_mpi_c 
        INTERFACE MPI::MPI_C)
    set_property(TARGET qe_mpi_c
        PROPERTY INTERFACE_COMPILE_DEFINITIONS ${mpi_definitions})
endif(QE_ENABLE_MPI)

###########################################################
# Tests
# Any executable target marked as test runner via
# 'add_test()' will be run by the 'test' make target
###########################################################
if(QE_ENABLE_TEST)
    enable_testing()
endif(QE_ENABLE_TEST)

###########################################################
# Components
###########################################################
add_subdirectory(external)
add_subdirectory(clib)
add_subdirectory(FFTXlib)
add_subdirectory(UtilXlib)
add_subdirectory(Modules)
add_subdirectory(LAXlib)
add_subdirectory(KS_Solvers)
add_subdirectory(dft-d3)
add_subdirectory(PW)
add_subdirectory(CPV)
add_subdirectory(atomic)

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

###########################################################
# Dependency graph generation
# Defines the custom target 'depgraph'
###########################################################
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/CMakeGraphVizOptions.cmake
               ${CMAKE_CURRENT_BINARY_DIR}/CMakeGraphVizOptions.cmake COPYONLY)
add_custom_target(depgraph
    "${CMAKE_COMMAND}" "--graphviz=depgraph" .
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
