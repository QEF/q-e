# Find the FFTW3 library
#
# Original version of this file:
#   Copyright (c) 2015, Wenzel Jakob
#   https://github.com/wjakob/layerlab/blob/master/cmake/FindFFTW.cmake, commit 4d58bfdc28891b4f9373dfe46239dda5a0b561c6
# Modifications:
#   Copyright (c) 2017, Patrick Bos
#
# Usage:
#   find_package(FFTW3 [REQUIRED] [QUIET] [COMPONENTS component1 ... componentX])
#
# It sets the following variables:
#   FFTW3_FOUND                  ... true if fftw is found on the system
#   FFTW3_[component]_FOUND  ... true if the component is found on the system (see components below)
#   FFTW3_LIBRARIES              ... full paths to all found fftw libraries
#   FFTW3_[component]_LIB        ... full path to one of the components (see below)
#   FFTW3_INCLUDE_DIRS           ... fftw include directory paths
#
# The following variables will be checked by the function
#   FFTW3_USE_STATIC_LIBS        ... if true, only static libraries are found, otherwise both static and shared.
#   FFTW3_ROOT                   ... if set, the libraries are exclusively searched under this path
#
# This package supports the following components:
#   FLOAT
#   DOUBLE
#   LONGDOUBLE
#   FLOAT_THREADS
#   DOUBLE_THREADS
#   LONGDOUBLE_THREADS
#   FLOAT_OPENMP
#   DOUBLE_OPENMP
#   LONGDOUBLE_OPENMP
#   FLOAT_MPI
#   DOUBLE_MPI

# Check whether to search static or dynamic libs
set(CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES})

if(${FFTW3_USE_STATIC_LIBS})
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV})
endif()

# Check common variables
set(FFTW3_DIRS $ENV{FFTW_HOME})
set(FFTW3_DIRS ${FFTW_DIRS} $ENV{FFTW_ROOT})
set(FFTW3_DIRS ${FFTW_DIRS} $ENV{FFTW_DIR})
set(FFTW3_DIRS ${FFTW_DIRS} $ENV{FFTW_PATH})
set(FFTW3_DIRS ${FFTW_DIRS} $ENV{FFTW3_ROOT})
set(FFTW3_DIRS ${FFTW_DIRS} $ENV{FFTW3_DIR})
set(FFTW3_DIRS ${FFTW_DIRS} $ENV{FFTW3_PATH})

if(FFTW3_DIRS)
  if("DOUBLE" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(
      FFTW3_DOUBLE
      NAMES "fftw3"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("DOUBLE_THREADS" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(
      FFTW3_DOUBLE_THREADS
      NAMES "fftw3_threads"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("DOUBLE_OPENMP" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(
      FFTW3_DOUBLE_OPENMP
      NAMES "fftw3_omp"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("DOUBLE_MPI" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(
      FFTW3_DOUBLE_MPI
      NAMES "fftw3_mpi"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("FLOAT" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(
      FFTW3_FLOAT
      NAMES "fftw3f"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("FLOAT_THREADS" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(
      FFTW3_FLOAT_THREADS
      NAMES "fftw3f_threads"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("FLOAT_OPENMP" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(
      FFTW3_FLOAT_OPENMP
      NAMES "fftw3f_omp"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("FLOAT_MPI" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(
      FFTW3_FLOAT_MPI
      NAMES "fftw3f_mpi"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("LONGDOUBLE" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(
      FFTW3_LONGDOUBLE
      NAMES "fftw3l"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("LONGDOUBLE_THREADS" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(
      FFTW3_LONGDOUBLE_THREADS
      NAMES "fftw3l_threads"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("LONGDOUBLE_OPENMP" IN_LIST FFTW3_FIND_COMPONENTS)
    find_library(FFTW3_LONGDOUBLE_OPENMP
      NAMES "fftw3l_omp"
      PATHS 
        ${FFTW3_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  find_path(FFTW3_INCLUDE_DIRS
    NAMES 
      "fftw3.h"
      "fftw3.f"
    PATHS 
      ${FFTW3_DIRS}
      $ENV{C_INCLUDE_PATH}
    PATH_SUFFIXES
      "include"
      "inc"
    NO_DEFAULT_PATH)
    
else()
  # Check if we can use PkgConfig
  find_package(PkgConfig)

  # Determine from PKG
  if(PKG_CONFIG_FOUND AND NOT FFTW3_ROOT)
    pkg_check_modules(PKG_FFTW3 QUIET "fftw3")
  endif()

  # Try to find in the LD_LIBRARY_PATH
  string(REPLACE ":" ";" LIBRARY_DIRS "$ENV{LD_LIBRARY_PATH}")

  if("DOUBLE" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(
    FFTW3_DOUBLE
    NAMES "fftw3"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("DOUBLE_THREADS" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(
    FFTW3_DOUBLE_THREADS
    NAMES "fftw3_threads"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("DOUBLE_OPENMP" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(
    FFTW3_DOUBLE_OPENMP
    NAMES "fftw3_omp"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS}
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("DOUBLE_MPI" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(
    FFTW3_DOUBLE_MPI
    NAMES "fftw3_mpi"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("FLOAT" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(
    FFTW3_FLOAT
    NAMES "fftw3f"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("FLOAT_THREADS" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(
    FFTW3_FLOAT_THREADS
    NAMES "fftw3f_threads"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS}
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("FLOAT_OPENMP" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(
    FFTW3_FLOAT_OPENMP
    NAMES "fftw3f_omp"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS}
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("FLOAT_MPI" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(
    FFTW3_FLOAT_MPI
    NAMES "fftw3f_mpi"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS}
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("LONGDOUBLE" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(
    FFTW3_LONGDOUBLE
    NAMES "fftw3l"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS}
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("LONGDOUBLE_THREADS" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(
    FFTW3_LONGDOUBLE_THREADS
    NAMES "fftw3l_threads"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("LONGDOUBLE_OPENMP" IN_LIST FFTW3_FIND_COMPONENTS)
  find_library(FFTW3_LONGDOUBLE_OPENMP
    NAMES "fftw3l_omp"
    PATHS 
      ${PKG_FFTW3_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  set(LD_DIRS)
  foreach(_lib ${LIBRARY_DIRS})
    get_filename_component(_root_lib_dir ${_lib} DIRECTORY)
    list(APPEND LD_DIRS ${_root_lib_dir})
  endforeach()
  
  find_path(FFTW3_INCLUDE_DIRS
    NAMES 
      "fftw3.h"
      "fftw3.f"
    PATHS 
      ${PKG_FFTW3_INCLUDE_DIRS} 
      ${INCLUDE_INSTALL_DIR}
      ${LD_DIRS}
      $ENV{C_INCLUDE_PATH}
    PATH_SUFFIXES
      "include"
      "inc")
endif()


# Components

add_library(FFTW3 INTERFACE IMPORTED)

if(FFTW3_DOUBLE)
  set(FFTW3_DOUBLE_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_DOUBLE})
else()
  set(FFTW3_DOUBLE_FOUND FALSE)
endif()

if(FFTW3_FLOAT)
  set(FFTW3_FLOAT_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_FLOAT})
else()
  set(FFTW3_FLOAT_FOUND FALSE)
endif()

if(FFTW3_LONGDOUBLE)
  set(FFTW3_LONGDOUBLE_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_LONGDOUBLE})
else()
  set(FFTW3_LONGDOUBLE_FOUND FALSE)
endif()

if(FFTW3_DOUBLE_THREADS)
  set(FFTW3_DOUBLE_THREADS_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_DOUBLE_THREADS})
else()
  set(FFTW3_DOUBLE_THREADS_FOUND FALSE)
endif()

if(FFTW3_FLOAT_THREADS)
  set(FFTW3_FLOAT_THREADS_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_FLOAT_THREADS})
else()
  set(FFTW3_FLOAT_THREADS_FOUND FALSE)
endif()

if(FFTW3_LONGDOUBLE_THREADS)
  set(FFTW3_LONGDOUBLE_THREADS_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_LONGDOUBLE_THREADS})
else()
  set(FFTW3_LONGDOUBLE_THREADS_FOUND FALSE)
endif()

if(FFTW3_DOUBLE_OPENMP)
  set(FFTW3_DOUBLE_OPENMP_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_DOUBLE_OPENMP})
  add_library(FFTW3::DoubleOpenMP INTERFACE IMPORTED)
else()
  set(FFTW3_DOUBLE_OPENMP_FOUND FALSE)
endif()

if(FFTW3_FLOAT_OPENMP)
  set(FFTW3_FLOAT_OPENMP_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_FLOAT_OPENMP})
else()
  set(FFTW3_FLOAT_OPENMP_FOUND FALSE)
endif()

if(FFTW3_LONGDOUBLE_OPENMP)
  set(FFTW3_LONGDOUBLE_OPENMP_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_LONGDOUBLE_OPENMP})
else()
  set(FFTW3_LONGDOUBLE_OPENMP_FOUND FALSE)
endif()

if(FFTW3_DOUBLE_MPI)
  set(FFTW3_DOUBLE_MPI_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_DOUBLE_MPI})
else()
  set(FFTW3_DOUBLE_MPI_FOUND FALSE)
endif()

if(FFTW3_FLOAT_MPI)
  set(FFTW3_FLOAT_MPI_FOUND TRUE)
  set(FFTW3_LIBRARIES ${FFTW3_LIBRARIES} ${FFTW3_FLOAT_MPI})
else()
  set(FFTW3_FLOAT_MPI_FOUND FALSE)
endif()

target_link_libraries(FFTW3 INTERFACE ${FFTW3_LIBRARIES})
target_include_directories(FFTW3 INTERFACE ${FFTW3_INCLUDE_DIRS})

# End components

set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3 REQUIRED_VARS FFTW3_LIBRARIES FFTW3_INCLUDE_DIRS)

mark_as_advanced(
        FFTW3_INCLUDE_DIRS
        FFTW3_LIBRARIES
        FFTW3_FLOAT
        FFTW3_DOUBLE
        FFTW3_LONGDOUBLE
        FFTW3_FLOAT_THREADS
        FFTW3_DOUBLE_THREADS
        FFTW3_LONGDOUBLE_THREADS
        FFTW3_FLOAT_OPENMP
        FFTW3_DOUBLE_OPENMP
        FFTW3_LONGDOUBLE_OPENMP
        FFTW3_FLOAT_MPI
        FFTW3_DOUBLE_MPI)