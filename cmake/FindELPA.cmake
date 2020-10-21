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
#  find_package(ELPA
#               [REQUIRED]             # Fail with error if elpa is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  ELPA depends on the following libraries:
#   - SCALAPACK
#
#  COMPONENTS are optional libraries ELPA could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - no components are available for now: maybe PLASMA in the future?
#
# Results are reported in variables:
#  ELPA_FOUND            - True if headers and requested libraries were found
#  ELPA_LINKER_FLAGS     - list of required linker flags (excluding -l and -L)
#  ELPA_INCLUDE_DIRS     - elpa include directories
#  ELPA_LIBRARY_DIRS     - Link directories for elpa libraries
#  ELPA_LIBRARIES        - elpa libraries
#  ELPA_INCLUDE_DIRS_DEP - elpa + dependencies include directories
#  ELPA_LIBRARY_DIRS_DEP - elpa + dependencies link directories
#  ELPA_LIBRARIES_DEP    - elpa libraries + dependencies
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DELPA_DIR=path/to/elpa):
#  ELPA_DIR              - Where to find the base directory of elpa
#  ELPA_INCDIR           - Where to find the header files
#  ELPA_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: ELPA_DIR, ELPA_INCDIR, ELPA_LIBDIR
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


if(NOT ELPA_FOUND)
  set(ELPA_DIR "" CACHE PATH "Installation directory of ELPA library")
  if (NOT ELPA_FIND_QUIETLY)
    message(STATUS "A cache variable, namely ELPA_DIR, has been set to specify the install directory of ELPA")
  endif()
endif(NOT ELPA_FOUND)

# ELPA depends on SCALAPACK anyway, try to find it
if (NOT SCALAPACK_FOUND)
  if(ELPA_FIND_REQUIRED)
    find_package(SCALAPACK REQUIRED)
  else()
    find_package(SCALAPACK)
  endif()
endif()


set(ENV_ELPA_DIR "$ENV{ELPA_DIR}")
set(ENV_ELPA_INCDIR "$ENV{ELPA_INCDIR}")
set(ENV_ELPA_LIBDIR "$ENV{ELPA_LIBDIR}")
set(ELPA_GIVEN_BY_USER "FALSE")
if ( ELPA_DIR OR ( ELPA_INCDIR AND ELPA_LIBDIR) OR ENV_ELPA_DIR OR (ENV_ELPA_INCDIR AND ENV_ELPA_LIBDIR) )
  set(ELPA_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT ELPA_GIVEN_BY_USER)

  pkg_search_module(ELPA elpa)
  if (NOT ELPA_FIND_QUIETLY)
    if (ELPA_FOUND AND ELPA_LIBRARIES)
      message(STATUS "Looking for ELPA - found using PkgConfig")
      #if(NOT ELPA_INCLUDE_DIRS)
      #    message("${Magenta}ELPA_INCLUDE_DIRS is empty using PkgConfig."
      #        "Perhaps the path to elpa headers is already present in your"
      #        "C(PLUS)_INCLUDE_PATH environment variable.${ColourReset}")
      #endif()
    else()
      message(STATUS "${Magenta}Looking for ELPA - not found using PkgConfig. "
	"\n   Perhaps you should add the directory containing elpa.pc "
	"\n   to the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (ELPA_FIND_VERSION_EXACT)
    if( NOT (ELPA_FIND_VERSION_MAJOR STREQUAL ELPA_VERSION_MAJOR) OR
	NOT (ELPA_FIND_VERSION_MINOR STREQUAL ELPA_VERSION_MINOR) )
      if(NOT ELPA_FIND_QUIETLY)
	message(FATAL_ERROR
	  "ELPA version found is ${ELPA_VERSION_STRING} "
	  "when required is ${ELPA_FIND_VERSION}")
      endif()
    endif()
  else()
    # if the version found is older than the required then error
    if( (ELPA_FIND_VERSION_MAJOR STRGREATER ELPA_VERSION_MAJOR) OR
	(ELPA_FIND_VERSION_MINOR STRGREATER ELPA_VERSION_MINOR) )
      if(NOT ELPA_FIND_QUIETLY)
	message(FATAL_ERROR
	  "ELPA version found is ${ELPA_VERSION_STRING} "
	  "when required is ${ELPA_FIND_VERSION} or newer")
      endif()
    endif()
  endif()

  # if pkg-config is used: these variables are empty
  # the pkg_search_module call will set the following:
  # ELPA_LDFLAGS: all required linker flags
  # ELPA_CFLAGS:  all required cflags
  set(ELPA_INCLUDE_DIRS_DEP "")
  set(ELPA_LIBRARY_DIRS_DEP "")
  set(ELPA_LIBRARIES_DEP "")
  # replace it anyway: we should update it with dependencies given by pkg-config
  set(ELPA_INCLUDE_DIRS_DEP "${ELPA_INCLUDE_DIRS}")
  set(ELPA_LIBRARY_DIRS_DEP "${ELPA_LIBRARY_DIRS}")
  set(ELPA_LIBRARIES_DEP "${ELPA_LIBRARIES}")

endif(PKG_CONFIG_EXECUTABLE AND NOT ELPA_GIVEN_BY_USER)

# if ELPA is not found using pkg-config
if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT ELPA_FOUND) OR (ELPA_GIVEN_BY_USER) )

  if (NOT ELPA_FIND_QUIETLY)
    message(STATUS "Looking for ELPA - PkgConfig not used")
  endif()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_ELPA_DIR "$ENV{ELPA_DIR}")
  set(ENV_ELPA_INCDIR "$ENV{ELPA_INCDIR}")
  if(ENV_ELPA_INCDIR)
    list(APPEND _inc_env "${ENV_ELPA_INCDIR}")
  elseif(ENV_ELPA_DIR)
    list(APPEND _inc_env "${ENV_ELPA_DIR}")
    list(APPEND _inc_env "${ENV_ELPA_DIR}/include")
    list(APPEND _inc_env "${ENV_ELPA_DIR}/include/elpa")
  else()
    if(WIN32)
      string(REPLACE ":" ";" _inc_env "$ENV{INCLUDE}")
    else()
      string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
      list(APPEND _inc_env "${_path_env}")
      string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
      list(APPEND _inc_env "${_path_env}")
      string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
      list(APPEND _inc_env "${_path_env}")
      string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
      list(APPEND _inc_env "${_path_env}")
    endif()
  endif()
  list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
  list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
  list(REMOVE_DUPLICATES _inc_env)


  # Try to find the elpa header in the given paths
  # -------------------------------------------------
  # call cmake macro to find the header path
  if(ELPA_INCDIR)
    set(ELPA_elpa.h_DIRS "ELPA_elpa.h_DIRS-NOTFOUND")
    find_path(ELPA_elpa.h_DIRS
      NAMES elpa.h
      HINTS ${ELPA_INCDIR})
  else()
    if(ELPA_DIR)
      set(ELPA_elpa.h_DIRS "ELPA_elpa.h_DIRS-NOTFOUND")
      find_path(ELPA_elpa.h_DIRS
	NAMES elpa.h
	HINTS ${ELPA_DIR}
	PATH_SUFFIXES "include" "include/elpa")
      if(NOT ELPA_elpa.h_DIRS)
        file(GLOB ELPA_elpa.h_PATH "${ELPA_DIR}/include/elpa-20*/elpa/elpa.h")
        get_filename_component(ELPA_elpa.h_DIRS "${ELPA_elpa.h_PATH}" PATH)
      endif()
    else()
      set(ELPA_elpa.h_DIRS "ELPA_elpa.h_DIRS-NOTFOUND")
      find_path(ELPA_elpa.h_DIRS
	NAMES elpa.h
	HINTS ${_inc_env})
    endif()
  endif()
  mark_as_advanced(ELPA_elpa.h_DIRS)

  # If found, add path to cmake variable
  # ------------------------------------
  if (ELPA_elpa.h_DIRS)
    message(STATUS "elpa.h found in ${ELPA_elpa.h_DIRS}")
    set(ELPA_INCLUDE_DIRS "${ELPA_elpa.h_DIRS}")
  else ()
    set(ELPA_INCLUDE_DIRS "ELPA_INCLUDE_DIRS-NOTFOUND")
    if(NOT ELPA_FIND_QUIETLY)
      message(STATUS "Looking for elpa -- elpa.h not found")
    endif()
  endif()

  # If not defined, guess the version string
  if(NOT ELPA_VERSION_STRING)
    string(REGEX MATCH "elpa-20[0-9][0-9]\.[0-9][0-9]\.[0-9][0-9][0-9]" CMAKE_MATCH_ELPA_VER "${ELPA_INCLUDE_DIRS}")
    string(REGEX MATCH "20[0-9][0-9]\.[0-9][0-9]\.[0-9][0-9][0-9]" ELPA_VERSION_STRING "${CMAKE_MATCH_ELPA_VER}")
    message(STATUS "ELPA version ${ELPA_VERSION_STRING}")
  endif()

  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_ELPA_LIBDIR "$ENV{ELPA_LIBDIR}")
  if(ENV_ELPA_LIBDIR)
    list(APPEND _lib_env "${ENV_ELPA_LIBDIR}")
  elseif(ENV_ELPA_DIR)
    list(APPEND _lib_env "${ENV_ELPA_DIR}")
    list(APPEND _lib_env "${ENV_ELPA_DIR}/lib")
  else()
    if(WIN32)
      string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
    else()
      if(APPLE)
	string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
      else()
	string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
      endif()
      list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
      list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
    endif()
  endif()
  list(REMOVE_DUPLICATES _lib_env)

  # Try to find the elpa lib in the given paths
  # ----------------------------------------------

  # call cmake macro to find the lib path
  if(ELPA_LIBDIR)
    set(ELPA_elpa_LIBRARY "ELPA_elpa_LIBRARY-NOTFOUND")
    find_library(ELPA_elpa_LIBRARY
      NAMES elpa
      HINTS ${ELPA_LIBDIR})
  else()
    if(ELPA_DIR)
      set(ELPA_elpa_LIBRARY "ELPA_elpa_LIBRARY-NOTFOUND")
      find_library(ELPA_elpa_LIBRARY
	NAMES elpa
	HINTS ${ELPA_DIR}
	PATH_SUFFIXES lib lib32 lib64)
    else()
      set(ELPA_elpa_LIBRARY "ELPA_elpa_LIBRARY-NOTFOUND")
      find_library(ELPA_elpa_LIBRARY
	NAMES elpa
	HINTS ${_lib_env})
    endif()
  endif()
  mark_as_advanced(ELPA_elpa_LIBRARY)

  # If found, add path to cmake variable
  # ------------------------------------
  if (ELPA_elpa_LIBRARY)
    get_filename_component(elpa_lib_path "${ELPA_elpa_LIBRARY}" PATH)
    # set cmake variables
    set(ELPA_LIBRARIES    "${ELPA_elpa_LIBRARY}")
    set(ELPA_LIBRARY_DIRS "${elpa_lib_path}")
  else ()
    set(ELPA_LIBRARIES    "ELPA_LIBRARIES-NOTFOUND")
    set(ELPA_LIBRARY_DIRS "ELPA_LIBRARY_DIRS-NOTFOUND")
    if(NOT ELPA_FIND_QUIETLY)
      message(STATUS "Looking for elpa -- lib elpa not found")
    endif()
  endif ()

  # check a function to validate the find
  if (ELPA_LIBRARIES)

    set(REQUIRED_LDFLAGS)
    set(REQUIRED_INCDIRS)
    set(REQUIRED_LIBDIRS)
    set(REQUIRED_LIBS)

    # ELPA
    if (ELPA_INCLUDE_DIRS)
      set(REQUIRED_INCDIRS "${ELPA_INCLUDE_DIRS}")
    endif()
    if (ELPA_LIBRARY_DIRS)
      set(REQUIRED_LIBDIRS "${ELPA_LIBRARY_DIRS}")
    endif()
    set(REQUIRED_LIBS "${ELPA_LIBRARIES}")
    # CBLAS
    if (CBLAS_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${CBLAS_INCLUDE_DIRS_DEP}")
    elseif (CBLAS_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${CBLAS_INCLUDE_DIRS}")
    endif()
    if(CBLAS_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${CBLAS_LIBRARY_DIRS_DEP}")
    elseif(CBLAS_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${CBLAS_LIBRARY_DIRS}")
    endif()
    if (CBLAS_LIBRARIES_DEP)
      list(APPEND REQUIRED_LIBS "${CBLAS_LIBRARIES_DEP}")
    elseif(CBLAS_LIBRARIES)
      list(APPEND REQUIRED_LIBS "${CBLAS_LIBRARIES}")
    endif()
    if (BLAS_LINKER_FLAGS)
      list(APPEND REQUIRED_LDFLAGS "${BLAS_LINKER_FLAGS}")
    endif()
    # LAPACK
    if (LAPACK_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${LAPACK_INCLUDE_DIRS}")
    endif()
    if(LAPACK_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${LAPACK_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${LAPACK_LIBRARIES}")
    if (LAPACK_LINKER_FLAGS)
      list(APPEND REQUIRED_LDFLAGS "${LAPACK_LINKER_FLAGS}")
    endif()
    # CUDA
    if (CUDA_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${CUDA_INCLUDE_DIRS}")
    endif()
    if(CUDA_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${CUDA_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${CUDA_CUBLAS_LIBRARIES};${CUDA_LIBRARIES}")

    # set required libraries for link
    set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
    set(CMAKE_REQUIRED_LIBRARIES)
    list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LDFLAGS}")
    foreach(lib_dir ${REQUIRED_LIBDIRS})
      list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
    endforeach()
    list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
    string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

    # test link
    unset(ELPA_WORKS CACHE)
    include(CheckFunctionExists)
    if(ELPA_VERSION_STRING AND ELPA_VERSION_STRING VERSION_LESS "2017")
      set(ELPA_TEST_FUNCTION get_elpa_row_col_comms)
    else()
      set(ELPA_TEST_FUNCTION elpa_init)
    endif()
    check_function_exists(${ELPA_TEST_FUNCTION} ELPA_WORKS)
    mark_as_advanced(ELPA_WORKS)

    if(ELPA_WORKS)
      # save link with dependencies
      set(ELPA_LIBRARIES_DEP    "${REQUIRED_LIBS}")
      set(ELPA_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
      set(ELPA_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
      set(ELPA_LINKER_FLAGS     "${REQUIRED_LDFLAGS}")
      list(REMOVE_DUPLICATES ELPA_LIBRARY_DIRS_DEP)
      list(REMOVE_DUPLICATES ELPA_INCLUDE_DIRS_DEP)
      list(REMOVE_DUPLICATES ELPA_LINKER_FLAGS)
    else()
      if(NOT ELPA_FIND_QUIETLY)
	message(STATUS "Looking for elpa : test of ${ELPA_TEST_FUNCTION} with
		elpa, cblas, cuda and lapack libraries fails")
	message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
	message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
	message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      endif()
    endif()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_FLAGS)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif(ELPA_LIBRARIES)

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT ELPA_FOUND) OR (ELPA_GIVEN_BY_USER) )

if (ELPA_LIBRARIES)
  if (ELPA_LIBRARY_DIRS)
    foreach(dir ${ELPA_LIBRARY_DIRS})
      if ("${dir}" MATCHES "elpa")
	set(first_lib_path "${dir}")
      endif()
    endforeach()
  else()
    list(GET ELPA_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(ELPA_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of ELPA library" FORCE)
  else()
    set(ELPA_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of ELPA library" FORCE)
  endif()
endif()
mark_as_advanced(ELPA_DIR)
mark_as_advanced(ELPA_DIR_FOUND)

# check that ELPA has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
if (PKG_CONFIG_EXECUTABLE AND ELPA_FOUND)
  find_package_handle_standard_args(ELPA DEFAULT_MSG
    ELPA_LIBRARIES)
else()
  find_package_handle_standard_args(ELPA DEFAULT_MSG
    ELPA_LIBRARIES
    ELPA_WORKS)
endif()
