###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find SCALAPACK library
# This module finds an installed fortran library that implements the SCALAPACK
# linear-algebra interface.
#
# This module sets the following variables:
#  SCALAPACK_FOUND - set to true if a library implementing the SCALAPACK interface
#    is found
#  SCALAPACK_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#    and -L).
#  SCALAPACK_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use SCALAPACK
#  SCALAPACK95_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use SCALAPACK95
#  SCALAPACK95_FOUND - set to true if a library implementing the SCALAPACK f95
#    interface is found
#  BLA_STATIC  if set on this determines what kind of linkage we do (static)
#  BLA_VENDOR  if set checks only the specified vendor, if not set checks
#     all the possibilities
#  BLA_F95     if set on tries to find the f95 interfaces for BLAS/SCALAPACK
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DSCALAPACK_DIR=path/to/scalapack):
#  SCALAPACK_DIR            - Where to find the base directory of scalapack
#  SCALAPACK_INCDIR         - Where to find the header files
#  SCALAPACK_LIBDIR         - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: SCALAPACK_DIR, SCALAPACK_INCDIR, SCALAPACK_LIBDIR
# Note that if BLAS_DIR is set, it will also look for scalapack in it
### List of vendors (BLA_VENDOR) valid in this module
##  Intel(mkl), ACML, Apple, NAS, Generic

#=============================================================================
# Copyright 2007-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)


## Some macros to print status when search for headers and libs
# This macro informs why the _lib_to_find file has not been found
macro(Print_Find_Library_Blas_Status _libname _lib_to_find)

  # save _libname upper/lower case
  string(TOUPPER ${_libname} LIBNAME)
  string(TOLOWER ${_libname} libname)

  # print status
  #message(" ")
  if(${LIBNAME}_LIBDIR)
    message("${Yellow}${LIBNAME}_LIBDIR is defined but ${_lib_to_find}"
      "has not been found in ${ARGN}${ColourReset}")
  else()
    if(${LIBNAME}_DIR)
      message("${Yellow}${LIBNAME}_DIR is defined but ${_lib_to_find}"
	"has not been found in ${ARGN}${ColourReset}")
    else()
      message("${Yellow}${_lib_to_find} not found."
	"Nor ${LIBNAME}_DIR neither ${LIBNAME}_LIBDIR"
	"are defined so that we look for ${_lib_to_find} in"
	"system paths (Linux: LD_LIBRARY_PATH, Windows: LIB,"
	"Mac: DYLD_LIBRARY_PATH,"
	"CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES,"
	"CMAKE_C_IMPLICIT_LINK_DIRECTORIES)${ColourReset}")
      if(_lib_env)
	message("${Yellow}${_lib_to_find} has not been found in"
	  "${_lib_env}${ColourReset}")
      endif()
    endif()
  endif()
  message("${BoldYellow}Please indicate where to find ${_lib_to_find}. You have three options:\n"
    "- Option 1: Provide the installation directory of the library with cmake option: -D${LIBNAME}_DIR=your/path/to/${libname}/\n"
    "- Option 2: Provide the directory where to find the library with cmake option: -D${LIBNAME}_LIBDIR=your/path/to/${libname}/lib/\n"
    "- Option 3: Update your environment variable (Linux: LD_LIBRARY_PATH, Windows: LIB, Mac: DYLD_LIBRARY_PATH)\n"
    "- Option 4: If your library provides a PkgConfig file, make sure pkg-config finds your library${ColourReset}")

endmacro()

if (NOT SCALAPACK_FOUND)
  set(SCALAPACK_DIR "" CACHE PATH "Installation directory of SCALAPACK library")
  if (NOT SCALAPACK_FIND_QUIETLY)
    message(STATUS "A cache variable, namely SCALAPACK_DIR, has been set to specify the install directory of SCALAPACK")
  endif()
endif (NOT SCALAPACK_FOUND)

option(SCALAPACK_VERBOSE "Print some additional information during SCALAPACK
libraries detection" OFF)
mark_as_advanced(SCALAPACK_VERBOSE)
if (BLAS_VERBOSE)
  set(SCALAPACK_VERBOSE ON)
endif ()
set(_scalapack_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
if (NOT _LANGUAGES_ MATCHES Fortran)
  include(CheckFunctionExists)
else (NOT _LANGUAGES_ MATCHES Fortran)
  include(CheckFortranFunctionExists)
endif (NOT _LANGUAGES_ MATCHES Fortran)

set(SCALAPACK_FOUND FALSE)
set(SCALAPACK95_FOUND FALSE)

# TODO: move this stuff to separate module

macro(Check_Scalapack_Libraries LIBRARIES _prefix _name _flags _list _blaslapack _mpi _threads)
  # This macro checks for the existence of the combination of fortran libraries
  # given by _list.  If the combination is found, this macro checks (using the
  # Check_Fortran_Function_Exists macro) whether can link against that library
  # combination using the name of a routine given by _name using the linker
  # flags given by _flags.  If the combination of libraries is found and passes
  # the link test, LIBRARIES is set to the list of complete library paths that
  # have been found.  Otherwise, LIBRARIES is set to FALSE.

  # N.B. _prefix is the prefix applied to the names of all cached variables that
  # are generated internally and marked advanced by this macro.

  set(_libraries_work TRUE)
  set(${LIBRARIES})
  set(_combined_name)
  set(ENV_MKLROOT "$ENV{MKLROOT}")
  set(ENV_BLAS_DIR "$ENV{BLAS_DIR}")
  set(ENV_BLAS_LIBDIR "$ENV{BLAS_LIBDIR}")
  set(ENV_SCALAPACK_DIR "$ENV{SCALAPACK_DIR}")
  set(ENV_SCALAPACK_LIBDIR "$ENV{SCALAPACK_LIBDIR}")
  if (NOT _libdir)
    if (BLAS_LIBDIR)
      list(APPEND _libdir "${BLAS_LIBDIR}")
    elseif (BLAS_DIR)
      list(APPEND _libdir "${BLAS_DIR}")
      list(APPEND _libdir "${BLAS_DIR}/lib")
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
	list(APPEND _libdir "${BLAS_DIR}/lib64")
	list(APPEND _libdir "${BLAS_DIR}/lib/intel64")
      else()
	list(APPEND _libdir "${BLAS_DIR}/lib32")
	list(APPEND _libdir "${BLAS_DIR}/lib/ia32")
      endif()
    elseif(ENV_BLAS_LIBDIR)
      list(APPEND _libdir "${ENV_BLAS_LIBDIR}")
    elseif(ENV_BLAS_DIR)
      list(APPEND _libdir "${ENV_BLAS_DIR}")
      list(APPEND _libdir "${ENV_BLAS_DIR}/lib")
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
	list(APPEND _libdir "${ENV_BLAS_DIR}/lib64")
	list(APPEND _libdir "${ENV_BLAS_DIR}/lib/intel64")
      else()
	list(APPEND _libdir "${ENV_BLAS_DIR}/lib32")
	list(APPEND _libdir "${ENV_BLAS_DIR}/lib/ia32")
      endif()
    endif()
    if (SCALAPACK_LIBDIR)
      list(APPEND _libdir "${SCALAPACK_LIBDIR}")
    elseif (SCALAPACK_DIR)
      list(APPEND _libdir "${SCALAPACK_DIR}")
      list(APPEND _libdir "${SCALAPACK_DIR}/lib")
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
	list(APPEND _libdir "${SCALAPACK_DIR}/lib64")
	list(APPEND _libdir "${SCALAPACK_DIR}/lib/intel64")
      else()
	list(APPEND _libdir "${SCALAPACK_DIR}/lib32")
	list(APPEND _libdir "${SCALAPACK_DIR}/lib/ia32")
      endif()
    elseif(ENV_SCALAPACK_LIBDIR)
      list(APPEND _libdir "${ENV_SCALAPACK_LIBDIR}")
    elseif(ENV_SCALAPACK_DIR)
      list(APPEND _libdir "${ENV_SCALAPACK_DIR}")
      list(APPEND _libdir "${ENV_SCALAPACK_DIR}/lib")
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
	list(APPEND _libdir "${ENV_SCALAPACK_DIR}/lib64")
	list(APPEND _libdir "${ENV_SCALAPACK_DIR}/lib/intel64")
      else()
	list(APPEND _libdir "${ENV_SCALAPACK_DIR}/lib32")
	list(APPEND _libdir "${ENV_SCALAPACK_DIR}/lib/ia32")
      endif()
    else()
      if (ENV_MKLROOT)
	list(APPEND _libdir "${ENV_MKLROOT}/lib")
	if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
	  list(APPEND _libdir "${ENV_MKLROOT}/lib64")
	  list(APPEND _libdir "${ENV_MKLROOT}/lib/intel64")
	else()
	  list(APPEND _libdir "${ENV_MKLROOT}/lib32")
	  list(APPEND _libdir "${ENV_MKLROOT}/lib/ia32")
	endif()
      endif()
      if (WIN32)
	string(REPLACE ":" ";" _libdir2 "$ENV{LIB}")
      elseif (APPLE)
	string(REPLACE ":" ";" _libdir2 "$ENV{DYLD_LIBRARY_PATH}")
      else ()
	string(REPLACE ":" ";" _libdir2 "$ENV{LD_LIBRARY_PATH}")
      endif ()
      list(APPEND _libdir "${_libdir2}")
      list(APPEND _libdir "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
      list(APPEND _libdir "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
    endif ()
  endif ()

  if (SCALAPACK_VERBOSE)
    message("${Cyan}Try to find SCALAPACK libraries: ${_list}")
  endif ()

  foreach(_library ${_list})
    set(_combined_name ${_combined_name}_${_library})

    if(_libraries_work)
      if (BLA_STATIC)
	if (WIN32)
	  set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
	endif ( WIN32 )
	if (APPLE)
	  set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
	else (APPLE)
	  set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
	endif (APPLE)
      else (BLA_STATIC)
	if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
	  # for ubuntu's libblas3gf and libscalapack3gf packages
	  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
	endif ()
      endif (BLA_STATIC)
      find_library(${_prefix}_${_library}_LIBRARY
	NAMES ${_library}
	HINTS ${_libdir}
	)
      mark_as_advanced(${_prefix}_${_library}_LIBRARY)
      # Print status if not found
      # -------------------------
      if (NOT ${_prefix}_${_library}_LIBRARY AND NOT SCALAPACK_FIND_QUIETLY AND SCALAPACK_VERBOSE)
	Print_Find_Library_Blas_Status(scalapack ${_library} ${_libdir})
      endif ()
      set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
      set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
    endif(_libraries_work)
  endforeach(_library ${_list})

  if(_libraries_work)
    # Test this combination of libraries.
    if(UNIX AND BLA_STATIC)
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} "-Wl,--start-group" ${${LIBRARIES}} ${_blaslapack} "-Wl,--end-group" ${_mpi} ${_threads})
    else(UNIX AND BLA_STATIC)
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_blaslapack} ${_mpi} ${_threads})
    endif(UNIX AND BLA_STATIC)
    if (SCALAPACK_VERBOSE)
      message("${Cyan}SCALAPACK libs found. Try to compile symbol ${_name} with"
	"following libraries: ${CMAKE_REQUIRED_LIBRARIES}")
    endif ()
    if(NOT SCALAPACK_FOUND)
      unset(${_prefix}${_combined_name}_WORKS CACHE)
    endif()
    if (NOT _LANGUAGES_ MATCHES Fortran)
      check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
    else (NOT _LANGUAGES_ MATCHES Fortran)
      check_fortran_function_exists(${_name} ${_prefix}${_combined_name}_WORKS)
    endif (NOT _LANGUAGES_ MATCHES Fortran)
    set(CMAKE_REQUIRED_LIBRARIES)
    mark_as_advanced(${_prefix}${_combined_name}_WORKS)
    set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  endif(_libraries_work)

  if(_libraries_work)
    set(${LIBRARIES} ${${LIBRARIES}} ${_blaslapack} ${_mpi} ${_threads})
  else(_libraries_work)
    set(${LIBRARIES} FALSE)
  endif(_libraries_work)

endmacro(Check_Scalapack_Libraries)


set(SCALAPACK_LINKER_FLAGS)
set(SCALAPACK_LIBRARIES)
set(SCALAPACK95_LIBRARIES)

if (NOT BLAS_FOUND)
  if(SCALAPACK_FIND_QUIETLY OR NOT SCALAPACK_FIND_REQUIRED)
    find_package(BLASEXT)
  else()
    find_package(BLASEXT REQUIRED)
  endif()
endif ()

if (NOT LAPACK_FOUND)
  if(SCALAPACK_FIND_QUIETLY OR NOT SCALAPACK_FIND_REQUIRED)
    find_package(LAPACKEXT)
  else()
    find_package(LAPACKEXT REQUIRED)
  endif()
endif ()

if (NOT MPI_FOUND)
  if(SCALAPACK_FIND_QUIETLY OR NOT SCALAPACK_FIND_REQUIRED)
    find_package(MPI)
  else()
    find_package(MPI REQUIRED)
  endif()
endif ()

if(BLAS_FOUND AND LAPACK_FOUND AND MPI_FOUND)
  set(SCALAPACK_LINKER_FLAGS ${BLAS_LINKER_FLAGS})
  list(APPEND SCALAPACK_LINKER_FLAGS ${LAPACK_LINKER_FLAGS})
  if ($ENV{BLA_VENDOR} MATCHES ".+")
    set(BLA_VENDOR $ENV{BLA_VENDOR})
  else ($ENV{BLA_VENDOR} MATCHES ".+")
    if(NOT BLA_VENDOR)
      set(BLA_VENDOR "All")
    endif(NOT BLA_VENDOR)
  endif ($ENV{BLA_VENDOR} MATCHES ".+")

  # Generic SCALAPACK library
  if (BLA_VENDOR STREQUAL "Generic" OR
      BLA_VENDOR STREQUAL "All")
    if ( NOT SCALAPACK_LIBRARIES )
      check_scalapack_libraries(
	SCALAPACK_LIBRARIES
	SCALAPACK
	pdgemm
	""
	"scalapack" # scalapack lib to look for
	"${LAPACK_LIBRARIES};${BLAS_LIBRARIES}" # blas and lapack libs
	"${MPI_Fortran_LIBRARIES}" # mpi libs
	""          # threads libs
	)
    endif ( NOT SCALAPACK_LIBRARIES )
  endif ()
  #intel scalapack
  if (BLA_VENDOR MATCHES "Intel" OR BLA_VENDOR STREQUAL "All")
    if (UNIX AND NOT WIN32)
      find_library(M_LIBRARY NAMES m)
      mark_as_advanced(M_LIBRARY)
      if(M_LIBRARY)
	set(LM "-lm")
      else()
	set(LM "")
      endif()
    endif ()
    if (_LANGUAGES_ MATCHES C OR _LANGUAGES_ MATCHES CXX)
      if(SCALAPACK_FIND_QUIETLY OR NOT SCALAPACK_FIND_REQUIRED)
	find_PACKAGE(Threads)
      else()
	find_package(Threads REQUIRED)
      endif()

      set(SCALAPACK_SEARCH_LIBS "")

      if (BLA_F95)
	set(SCALAPACK_mkl_SEARCH_SYMBOL "PDGEMM")
	set(_LIBRARIES SCALAPACK95_LIBRARIES)
	set(_BLAS_LIBRARIES ${BLAS95_LIBRARIES})
	list(APPEND SCALAPACK_SEARCH_LIBS "mkl_scalapack_lp64")
      else()
	set(SCALAPACK_mkl_SEARCH_SYMBOL "pdgemm")
	set(_LIBRARIES SCALAPACK_LIBRARIES)
	set(_BLAS_LIBRARIES ${BLAS_LIBRARIES})
	list(APPEND SCALAPACK_SEARCH_LIBS "mkl_scalapack_lp64")
      endif()

      # First try empty scalapack libs
      if (NOT ${_LIBRARIES})
	check_scalapack_libraries(
	  ${_LIBRARIES}
	  BLAS
	  ${SCALAPACK_mkl_SEARCH_SYMBOL}
	  ""
	  ""
	  "${_BLAS_LIBRARIES}"
	  ""
	  "${MPI_Fortran_LIBRARIES}"
	  )
      endif ()
      # Then try the search libs
      foreach (IT ${SCALAPACK_SEARCH_LIBS})
	if (NOT ${_LIBRARIES})
	  check_scalapack_libraries(
	    ${_LIBRARIES}
	    BLAS
	    ${SCALAPACK_mkl_SEARCH_SYMBOL}
	    ""
	    "${IT};mkl_blacs_intelmpi_lp64"
	    "${_BLAS_LIBRARIES}"
	    ""
	    "${MPI_Fortran_LIBRARIES}"
	    )
	endif ()
      endforeach ()
    endif ()
  endif()
else(BLAS_FOUND AND LAPACK_FOUND AND MPI_FOUND)
  message(STATUS "SCALAPACK requires BLAS, LAPACK, and MPI")
endif(BLAS_FOUND AND LAPACK_FOUND AND MPI_FOUND)

if(BLA_F95)
  if(SCALAPACK95_LIBRARIES)
    set(SCALAPACK95_FOUND TRUE)
  else(SCALAPACK95_LIBRARIES)
    set(SCALAPACK95_FOUND FALSE)
  endif(SCALAPACK95_LIBRARIES)
  if(NOT SCALAPACK_FIND_QUIETLY)
    if(SCALAPACK95_FOUND)
      message(STATUS "A library with SCALAPACK95 API found.")
      message(STATUS "SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARIES}")
    else(SCALAPACK95_FOUND)
      message(WARNING "BLA_VENDOR has been set to ${BLA_VENDOR} but SCALAPACK 95 libraries could not be found or check of symbols failed."
	"\nPlease indicate where to find SCALAPACK libraries. You have three options:\n"
	"- Option 1: Provide the installation directory of SCALAPACK library with cmake option: -DSCALAPACK_DIR=your/path/to/scalapack\n"
	"- Option 2: Provide the directory where to find BLAS libraries with cmake option: -DBLAS_LIBDIR=your/path/to/blas/libs\n"
	"- Option 3: Update your environment variable (Linux: LD_LIBRARY_PATH, Windows: LIB, Mac: DYLD_LIBRARY_PATH)\n"
	"\nTo follow libraries detection more precisely you can activate a verbose mode with -DSCALAPACK_VERBOSE=ON at cmake configure."
	"\nYou could also specify a BLAS vendor to look for by setting -DBLA_VENDOR=blas_vendor_name."
	"\nList of possible BLAS vendor: Goto, ATLAS PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_32 (intel mkl v10 32 bit),"
	"Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model), Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),"
	"Intel( older versions of mkl 32 and 64 bit), ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")
      if(SCALAPACK_FIND_REQUIRED)
	message(FATAL_ERROR
	  "A required library with SCALAPACK95 API not found. Please specify library location."
	  )
      else(SCALAPACK_FIND_REQUIRED)
	message(STATUS
	  "A library with SCALAPACK95 API not found. Please specify library location."
	  )
      endif(SCALAPACK_FIND_REQUIRED)
    endif(SCALAPACK95_FOUND)
  endif(NOT SCALAPACK_FIND_QUIETLY)
  set(SCALAPACK_FOUND "${SCALAPACK95_FOUND}")
  set(SCALAPACK_LIBRARIES "${SCALAPACK95_LIBRARIES}")
else(BLA_F95)
  if(SCALAPACK_LIBRARIES)
    set(SCALAPACK_FOUND TRUE)
  else(SCALAPACK_LIBRARIES)
    set(SCALAPACK_FOUND FALSE)
  endif(SCALAPACK_LIBRARIES)

  if(NOT SCALAPACK_FIND_QUIETLY)
    if(SCALAPACK_FOUND)
      message(STATUS "A library with SCALAPACK API found.")
      message(STATUS "SCALAPACK_LIBRARIES ${SCALAPACK_LIBRARIES}")
    else(SCALAPACK_FOUND)
      message(WARNING "BLA_VENDOR has been set to ${BLA_VENDOR} but SCALAPACK libraries could not be found or check of symbols failed."
	"\nPlease indicate where to find SCALAPACK libraries. You have three options:\n"
	"- Option 1: Provide the installation directory of SCALAPACK library with cmake option: -DSCALAPACK_DIR=your/path/to/scalapack\n"
	"- Option 2: Provide the directory where to find BLAS libraries with cmake option: -DBLAS_LIBDIR=your/path/to/blas/libs\n"
	"- Option 3: Update your environment variable (Linux: LD_LIBRARY_PATH, Windows: LIB, Mac: DYLD_LIBRARY_PATH)\n"
	"\nTo follow libraries detection more precisely you can activate a verbose mode with -DSCALAPACK_VERBOSE=ON at cmake configure."
	"\nYou could also specify a BLAS vendor to look for by setting -DBLA_VENDOR=blas_vendor_name."
	"\nList of possible BLAS vendor: Goto, ATLAS PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_32 (intel mkl v10 32 bit),"
	"Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model), Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),"
	"Intel( older versions of mkl 32 and 64 bit), ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")
      if(SCALAPACK_FIND_REQUIRED)
	message(FATAL_ERROR
	  "A required library with SCALAPACK API not found. Please specify library location."
	  )
      else(SCALAPACK_FIND_REQUIRED)
	message(STATUS
	  "A library with SCALAPACK API not found. Please specify library location."
	  )
      endif(SCALAPACK_FIND_REQUIRED)
    endif(SCALAPACK_FOUND)
  endif(NOT SCALAPACK_FIND_QUIETLY)
endif(BLA_F95)

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_scalapack_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

if (SCALAPACK_LIBRARIES)
  list(GET SCALAPACK_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (${first_lib_path} MATCHES "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)")
    string(REGEX REPLACE "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)" "" not_cached_dir "${first_lib_path}")
    set(SCALAPACK_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of SCALAPACK library" FORCE)
  else()
    set(SCALAPACK_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of SCALAPACK library" FORCE)
  endif()
endif()
mark_as_advanced(SCALAPACK_DIR)
mark_as_advanced(SCALAPACK_DIR_FOUND)
