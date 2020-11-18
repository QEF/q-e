# Find the FFTW library
#
# Original version of this file:
#   Copyright (c) 2015, Wenzel Jakob
#   https://github.com/wjakob/layerlab/blob/master/cmake/FindFFTW.cmake, commit 4d58bfdc28891b4f9373dfe46239dda5a0b561c6
# Modifications:
#   Copyright (c) 2017, Patrick Bos
#
# Usage:
#   find_package(FFTW [REQUIRED] [QUIET] [COMPONENTS component1 ... componentX])
#
# It sets the following variables:
#   FFTW_FOUND                  ... true if fftw is found on the system
#   FFTW_[component]_FOUND  ... true if the component is found on the system (see components below)
#   FFTW_LIBRARIES              ... full paths to all found fftw libraries
#   FFTW_[component]_LIB        ... full path to one of the components (see below)
#   FFTW_INCLUDE_DIRS           ... fftw include directory paths
#
# The following variables will be checked by the function
#   FFTW_USE_STATIC_LIBS        ... if true, only static libraries are found, otherwise both static and shared.
#   FFTW_ROOT                   ... if set, the libraries are exclusively searched under this path
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

if(${FFTW_USE_STATIC_LIBS})
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
else()
  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV})
endif()

# Check common variables
set(FFTW_DIRS $ENV{FFTW_HOME})
set(FFTW_DIRS ${FFTW_DIRS} $ENV{FFTW_ROOT})
set(FFTW_DIRS ${FFTW_DIRS} $ENV{FFTW_DIR})
set(FFTW_DIRS ${FFTW_DIRS} $ENV{FFTW_PATH})

if(FFTW_DIRS)
  if("DOUBLE" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(
      FFTW_DOUBLE
      NAMES "fftw3"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("DOUBLE_THREADS" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(
      FFTW_DOUBLE_THREADS
      NAMES "fftw3_threads"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("DOUBLE_OPENMP" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(
      FFTW_DOUBLE_OPENMP
      NAMES "fftw3_omp"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("DOUBLE_MPI" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(
      FFTW_DOUBLE_MPI
      NAMES "fftw3_mpi"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("FLOAT" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(
      FFTW_FLOAT
      NAMES "fftw3f"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("FLOAT_THREADS" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(
      FFTW_FLOAT_THREADS
      NAMES "fftw3f_threads"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("FLOAT_OPENMP" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(
      FFTW_FLOAT_OPENMP
      NAMES "fftw3f_omp"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("FLOAT_MPI" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(
      FFTW_FLOAT_MPI
      NAMES "fftw3f_mpi"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("LONGDOUBLE" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(
      FFTW_LONGDOUBLE
      NAMES "fftw3l"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("LONGDOUBLE_THREADS" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(
      FFTW_LONGDOUBLE_THREADS
      NAMES "fftw3l_threads"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  if("LONGDOUBLE_OPENMP" IN_LIST FFTW_FIND_COMPONENTS)
    find_library(FFTW_LONGDOUBLE_OPENMP
      NAMES "fftw3l_omp"
      PATHS 
        ${FFTW_DIRS}
      PATH_SUFFIXES
        "lib"
        "lib64"
      NO_DEFAULT_PATH)
  endif()

  find_path(FFTW_INCLUDE_DIRS
    NAMES 
      "fftw3.h"
      "fftw3.f"
    PATHS 
      ${FFTW_DIRS}
      $ENV{C_INCLUDE_PATH}
    PATH_SUFFIXES
      "include"
      "inc"
    NO_DEFAULT_PATH)
    
else()
  # Check if we can use PkgConfig
  find_package(PkgConfig)

  # Determine from PKG
  if(PKG_CONFIG_FOUND AND NOT FFTW_ROOT)
    pkg_check_modules(PKG_FFTW QUIET "fftw3")
  endif()

  # Try to find in the LD_LIBRARY_PATH
  string(REPLACE ":" ";" LIBRARY_DIRS "$ENV{LD_LIBRARY_PATH}")

  if("DOUBLE" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(
    FFTW_DOUBLE
    NAMES "fftw3"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("DOUBLE_THREADS" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(
    FFTW_DOUBLE_THREADS
    NAMES "fftw3_threads"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("DOUBLE_OPENMP" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(
    FFTW_DOUBLE_OPENMP
    NAMES "fftw3_omp"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS}
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("DOUBLE_MPI" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(
    FFTW_DOUBLE_MPI
    NAMES "fftw3_mpi"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("FLOAT" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(
    FFTW_FLOAT
    NAMES "fftw3f"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("FLOAT_THREADS" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(
    FFTW_FLOAT_THREADS
    NAMES "fftw3f_threads"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS}
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("FLOAT_OPENMP" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(
    FFTW_FLOAT_OPENMP
    NAMES "fftw3f_omp"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS}
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("FLOAT_MPI" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(
    FFTW_FLOAT_MPI
    NAMES "fftw3f_mpi"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS}
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("LONGDOUBLE" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(
    FFTW_LONGDOUBLE
    NAMES "fftw3l"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS}
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("LONGDOUBLE_THREADS" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(
    FFTW_LONGDOUBLE_THREADS
    NAMES "fftw3l_threads"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS} 
      ${LIB_INSTALL_DIR}
      ${LIBRARY_DIRS}
    PATH_SUFFIXES
      "lib"
      "lib64")
  endif()

  if("LONGDOUBLE_OPENMP" IN_LIST FFTW_FIND_COMPONENTS)
  find_library(FFTW_LONGDOUBLE_OPENMP
    NAMES "fftw3l_omp"
    PATHS 
      ${PKG_FFTW_LIBRARY_DIRS} 
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
  
  find_path(FFTW_INCLUDE_DIRS
    NAMES 
      "fftw3.h"
      "fftw3.f"
    PATHS 
      ${PKG_FFTW_INCLUDE_DIRS} 
      ${INCLUDE_INSTALL_DIR}
      ${LD_DIRS}
      $ENV{C_INCLUDE_PATH}
    PATH_SUFFIXES
      "include"
      "inc")
endif()


# Components

add_library(FFTW INTERFACE IMPORTED)

if(FFTW_DOUBLE)
  set(FFTW_DOUBLE_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE})
else()
  set(FFTW_DOUBLE_FOUND FALSE)
endif()

if(FFTW_FLOAT)
  set(FFTW_FLOAT_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_FLOAT})
else()
  set(FFTW_FLOAT_FOUND FALSE)
endif()

if(FFTW_LONGDOUBLE)
  set(FFTW_LONGDOUBLE_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_LONGDOUBLE})
else()
  set(FFTW_LONGDOUBLE_FOUND FALSE)
endif()

if(FFTW_DOUBLE_THREADS)
  set(FFTW_DOUBLE_THREADS_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_THREADS})
else()
  set(FFTW_DOUBLE_THREADS_FOUND FALSE)
endif()

if(FFTW_FLOAT_THREADS)
  set(FFTW_FLOAT_THREADS_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_FLOAT_THREADS})
else()
  set(FFTW_FLOAT_THREADS_FOUND FALSE)
endif()

if(FFTW_LONGDOUBLE_THREADS)
  set(FFTW_LONGDOUBLE_THREADS_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_LONGDOUBLE_THREADS})
else()
  set(FFTW_LONGDOUBLE_THREADS_FOUND FALSE)
endif()

if(FFTW_DOUBLE_OPENMP)
  set(FFTW_DOUBLE_OPENMP_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_OPENMP})
  add_library(FFTW::DoubleOpenMP INTERFACE IMPORTED)
else()
  set(FFTW_DOUBLE_OPENMP_FOUND FALSE)
endif()

if(FFTW_FLOAT_OPENMP)
  set(FFTW_FLOAT_OPENMP_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_FLOAT_OPENMP})
else()
  set(FFTW_FLOAT_OPENMP_FOUND FALSE)
endif()

if(FFTW_LONGDOUBLE_OPENMP)
  set(FFTW_LONGDOUBLE_OPENMP_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_LONGDOUBLE_OPENMP})
else()
  set(FFTW_LONGDOUBLE_OPENMP_FOUND FALSE)
endif()

if(FFTW_DOUBLE_MPI)
  set(FFTW_DOUBLE_MPI_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_DOUBLE_MPI})
else()
  set(FFTW_DOUBLE_MPI_FOUND FALSE)
endif()

if(FFTW_FLOAT_MPI)
  set(FFTW_FLOAT_MPI_FOUND TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARIES} ${FFTW_FLOAT_MPI})
else()
  set(FFTW_FLOAT_MPI_FOUND FALSE)
endif()

target_link_libraries(FFTW INTERFACE ${FFTW_LIBRARIES})
target_include_directories(FFTW INTERFACE ${FFTW_INCLUDE_DIRS})

# End components

set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW REQUIRED_VARS FFTW_LIBRARIES FFTW_INCLUDE_DIRS)

mark_as_advanced(
        FFTW_INCLUDE_DIRS
        FFTW_LIBRARIES
        FFTW_FLOAT
        FFTW_DOUBLE
        FFTW_LONGDOUBLE
        FFTW_FLOAT_THREADS
        FFTW_DOUBLE_THREADS
        FFTW_LONGDOUBLE_THREADS
        FFTW_FLOAT_OPENMP
        FFTW_DOUBLE_OPENMP
        FFTW_LONGDOUBLE_OPENMP
        FFTW_FLOAT_MPI
        FFTW_DOUBLE_MPI)