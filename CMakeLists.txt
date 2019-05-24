cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(qe
    VERSION 6.4.1
    LANGUAGES Fortran C
    HOMEPAGE_URL "https://www.quantum-espresso.org/"
    DESCRIPTION "ESPRESSO: opEn-Source Package for Research in Electronic Structure, Simulation, and Optimization")

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

###########################################################
# Build Options
###########################################################

option(QE_ENABLE_MPI
    "enable distributed execution support via MPI" ON)

# Each one of the following options boils down to a simple
# compile definition, no additional dependencies are needed
#
# TODO move all feature definitions to a proper config.h
option(QE_ENABLE_TRACE
    "enable execution tracing output" OFF)
option(QE_ENABLE_MPI_INPLACE
    "enable inplace MPI calls (ignored when QE_ENABLE_MPI=OFF)" OFF)
option(QE_ENABLE_BARRIER
    "enable global synchronization between execution units" OFF)

###########################################################
# Build Features
# The following targets will be defined:
add_library(qe_features INTERFACE)
add_library(QE::Features ALIAS qe_features)
###########################################################
# TODO change all ifdefs throughout code base to match
# cmake options
if(QE_ENABLE_TRACE)
    set_property(TARGET qe_features
        PROPERTY
            INTERFACE_COMPILE_DEFINITIONS __TRACE)
endif()
if(QE_ENABLE_MPI_INPLACE)
    set_property(TARGET qe_features
        PROPERTY
            INTERFACE_COMPILE_DEFINITIONS __USE_INPLACE_MPI)
endif()
if(QE_ENABLE_BARRIER)
    set_property(TARGET qe_features
        PROPERTY
            INTERFACE_COMPILE_DEFINITIONS __USE_BARRIER)
endif()

###########################################################
# MPI
# The following targets will be defined:
add_library(qe_mpi_fortran INTERFACE)
add_library(QE::MPI_Fortran ALIAS qe_mpi_fortran)
add_library(qe_mpi_c INTERFACE)
add_library(QE::MPI_C ALIAS qe_mpi_c)
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
# Utility functions
###########################################################

function(add_fortran_library LIB)
    add_library(${LIB} ${ARGN})
    # set module path to LIB_DIR/mod
    get_target_property(LIB_DIR ${LIB} BINARY_DIR)
    set_target_properties(${LIB} PROPERTIES Fortran_MODULE_DIRECTORY ${LIB_DIR}/mod)
    # making LIB_DIR/mod available for libraries linking LIB 
    target_include_directories(${LIB} INTERFACE ${LIB_DIR}/mod)
endfunction(add_fortran_library)

###########################################################
# Components
###########################################################

add_subdirectory(UtilXLib)
add_subdirectory(LAXLib)
