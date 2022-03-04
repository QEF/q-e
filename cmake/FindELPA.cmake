###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
# @copyright (c) 2017 King Abdullah University of Science and Technology. All rights reserved.
#
###
#
# - Find ELPA include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(ELPA [version] [EXACT] [QUIET] [REQUIRED] [COMPONENTS <comp1> ..])
#
#  ELPA depends on the following libraries:
#   - SCALAPACK
#
#  COMPONENTS are optional libraries ELPA could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - no components are available for now: maybe OpenMP CUDA in the future?
#
# Results are reported in variables:
#  ELPA_FOUND            - True if headers and requested libraries were found
#  ELPA_VERSION          - Version string
#  ELPA_INCLUDE_DIRS     - elpa include directories
#  ELPA_Fortran_MODS_DIR - elpa directory containing Fortran modules
#  ELPA_LIBRARIES        - elpa libraries
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DELPA_ROOT=path/to/elpa):
#  ELPA_ROOT             - Where to find the base directory of elpa
#  ELPA_DIR(deprecated)  - Where to find the base directory of elpa
# The module can also look for the following environment variables if paths
# are not given as cmake variable: ELPA_ROOT, ELPA_DIR(deprecated)
#
# searching precedence
# ELPA_ROOT > ENV{ELPA_ROOT} > CMAKE_INCLUDE_PATH/CMAKE_LIBRARY_PATH > CMAKE_PREFIX_PATH
# > ENV{CMAKE_INCLUDE_PATH}/ENV{CMAKE_LIBRARY_PATH}, ENV{CMAKE_PREFIX_PATH}
# > PKG_ELPA_INCLUDE_DIRS/PKG_ELPA_LIBRARY_DIRS > system locations
#
#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013-2016 Florent Pruvost
# Copyright 2017      Eduardo Gonzalez
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file ECRC-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Ecrc, substitute the full
#  License text for the above reference.)

# Check whether to search static or dynamic libs
set(CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES})

if(${ELPA_USE_STATIC_LIBS})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV})
endif()

# ELPA depends on SCALAPACK anyway, try to find it
if(NOT SCALAPACK_FOUND)
    if(ELPA_FIND_REQUIRED)
        find_package(SCALAPACK REQUIRED)
    else()
        find_package(SCALAPACK)
    endif()
endif()

# minimal support of ELPA_DIR ENV{ELPA_DIR}
if(NOT ELPA_ROOT AND NOT ENV{ELPA_ROOT} AND ELPA_DIR)
    message(WARNING "CMake variable ELPA_DIR has been deprecated")
    set(ELPA_ROOT ${ELPA_DIR})
endif()
if(NOT ELPA_ROOT AND NOT ENV{ELPA_ROOT} AND ENV{ELPA_DIR})
    message(WARNING "Environment variable ELPA_DIR has been deprecated")
    set(ELPA_ROOT $ENV{ELPA_DIR})
endif()

# When ELPA is manually installed, the header C file elpa.h is placed
# under the directory INSTALL_PREFIX/include/elpa-20XX/elpa
# In C source code, an include line uses "elpa/elpa.h".
# The include directory should be INSTALL_PREFIX/include/elpa-20XX
# glob_elpa_header_file handles versioned include directories
macro(glob_elpa_header_file)
    string(REPLACE ":" ";" ARGN_EXPANDED "${ARGN}")
    foreach(DIR ${ARGN_EXPANDED})
        message(DEBUG "globing ${DIR}")
        if(NOT ELPA_elpa.h_FILEPATH AND DIR)
            file(GLOB ELPA_elpa.h_FILEPATH "${DIR}/include/elpa-20*/elpa/elpa.h")
        endif()
    endforeach()
endmacro()

glob_elpa_header_file("${ELPA_ROOT}" "$ENV{ELPA_ROOT}" "${CMAKE_INCLUDE_PATH}" "${CMAKE_PREFIX_PATH}"
                      "$ENV{CMAKE_INCLUDE_PATH}" "$ENV{CMAKE_PREFIX_PATH}")

if(ELPA_elpa.h_FILEPATH)
    string(REGEX REPLACE "/elpa/elpa.h$" "" ELPA_elpa.h_INCLUDE_DIR "${ELPA_elpa.h_FILEPATH}")
    if(NOT ELPA_FIND_QUIETLY)
        message(STATUS "Globing elpa.h on versioned paths found ${ELPA_elpa.h_FILEPATH}")
    endif()
    find_path(
        ELPA_INCLUDE_DIRS
        NAMES "elpa/elpa.h"
        HINTS ${ELPA_elpa.h_INCLUDE_DIR}
        NO_DEFAULT_PATH)
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
if(NOT ELPA_INCLUDE_DIRS)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        if(NOT PKG_ELPA_FOUND AND ELPA_PKGCONFIG_VERSION)
            pkg_search_module(PKG_ELPA elpa-${ELPA_PKGCONFIG_VERSION})
        endif()
        if(NOT PKG_ELPA_FOUND)
            pkg_search_module(PKG_ELPA elpa)
        endif()
        if(NOT ELPA_FIND_QUIETLY)
            if(PKG_ELPA_FOUND)
                message(STATUS "Looking for ELPA - info found by PkgConfig")
                message("     PkgConfig found include directory: ${PKG_ELPA_INCLUDE_DIRS}")
                message("     PkgConfig found library directory: ${PKG_ELPA_LIBRARY_DIRS}")
            else()
                message(STATUS "Looking for ELPA - info not found by PkgConfig")
            endif()
        endif()
    endif()
endif()

# search for "elpa/elpa.h" again
# ELPA_INCLUDE_DIRS is a cached variable, if it has been found before, this find_path is no-op
find_path(
    ELPA_INCLUDE_DIRS
    NAMES "elpa/elpa.h"
    HINTS ${PKG_ELPA_INCLUDE_DIRS}
    PATH_SUFFIXES "include" "include/elpa")

if(ELPA_INCLUDE_DIRS)
    # Looking for module directory based on elpa.mod
    find_path(
        ELPA_Fortran_MODS_DIR
        NAMES elpa.mod
        HINTS ${ELPA_INCLUDE_DIRS}
        PATH_SUFFIXES "modules" NO_DEFAULT_PATH)

    # Looking for elpa1.mod if elpa.mod was not found
    find_path(
        ELPA_Fortran_MODS_DIR
        NAMES elpa1.mod
        HINTS ${ELPA_INCLUDE_DIRS}
        PATH_SUFFIXES "modules" NO_DEFAULT_PATH)
endif()

find_library(
    ELPA_LIBRARIES
    NAMES elpa
    HINTS ${PKG_ELPA_LIBRARY_DIRS}
    PATH_SUFFIXES "lib" "lib64")

# extract version string from ELPA_INCLUDE_DIRS
if(ELPA_INCLUDE_DIRS)
    string(REGEX MATCH "elpa-20[0-9][0-9]\.[0-9][0-9]\.[0-9][0-9][0-9]" CMAKE_MATCH_ELPA_VER "${ELPA_INCLUDE_DIRS}")
    string(REGEX MATCH "20[0-9][0-9]\.[0-9][0-9]\.[0-9][0-9][0-9]" ELPA_VERSION_FROM_INCLUDE_DIRS
                 "${CMAKE_MATCH_ELPA_VER}")
endif()
if(NOT ELPA_FIND_QUIETLY)
    if(ELPA_VERSION_FROM_INCLUDE_DIRS)
        message(STATUS "ELPA version string extracted from ELPA_INCLUDE_DIRS : ${ELPA_VERSION_FROM_INCLUDE_DIRS}")
    else()
        message(STATUS "Extracting ELPA version string from ELPA_INCLUDE_DIRS failed.")
    endif()
endif()

# If ELPA_VERSION is not defined, use ELPA_VERSION_FROM_INCLUDE_DIRS
if(ELPA_VERSION_FROM_INCLUDE_DIRS)
    if(NOT ELPA_VERSION)
        message(DEBUG "setting ELPA_VERSION by ELPA_VERSION_FROM_INCLUDE_DIRS")
        set(ELPA_VERSION ${ELPA_VERSION_FROM_INCLUDE_DIRS})
    elseif(NOT ELPA_VERSION_FROM_INCLUDE_DIRS STREQUAL ELPA_VERSION)
        message(
            WARNING "Confusing elpa versions exist in the environement. "
                    "Found ${ELPA_VERSION_FROM_INCLUDE_DIRS} through ELPA_INCLUDE_DIRS "
                    "but version string ${ELPA_VERSION} is in use. "
                    "Check carefully if the desired version got picked up.")
    endif()
endif()

# If ELPA_VERSION is not defined, use PKG_ELPA_VERSION
if(PKG_ELPA_VERSION)
    if(NOT ELPA_VERSION)
        message(DEBUG "setting ELPA_VERSION by PKG_ELPA_VERSION")
        set(ELPA_VERSION ${PKG_ELPA_VERSION})
    elseif(NOT PKG_ELPA_VERSION STREQUAL ELPA_VERSION)
        message(WARNING "Confusing elpa versions exist in the environement. "
                        "pkg-config found ${PKG_ELPA_VERSION} but version string ${ELPA_VERSION} is in use. "
                        "Check carefully if the desired version got picked up.")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    ELPA
    REQUIRED_VARS ELPA_LIBRARIES ELPA_INCLUDE_DIRS ELPA_Fortran_MODS_DIR ELPA_VERSION
    VERSION_VAR ELPA_VERSION)
