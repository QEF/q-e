#[=======================================================================[.rst:
FindVENDOR
--------

Find Vendor library

Input Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``VENDOR_STATIC``
  if ``ON`` use static linkage

``VENDOR``
  If set, checks only the specified vendor, if not set checks all the
  possibilities.  List of vendors valid in this module:

  * ``Intel10_32``        (intel mkl v10 32 bit)
  * ``Intel10_64lp``      (intel mkl v10+ 64 bit, threaded code, lp64 model)
  * ``Intel10_64lp_seq``  (intel mkl v10+ 64 bit, sequential code, lp64 model)
  * ``Intel10_64ilp``     (intel mkl v10+ 64 bit, threaded code, ilp64 model)
  * ``Intel10_64ilp_seq`` (intel mkl v10+ 64 bit, sequential code, ilp64 model)
  * ``Intel10_64_dyn``    (intel mkl v10+ 64 bit, single dynamic library)
  * ``Intel``             (obsolete versions of mkl 32 and 64 bit)
  * ``Armpl_64lp_mp``     (armpl, threaded code, lp64 model)
  * ``Armpl_64lp``        (armpl, sequential code, lp64 model)
  * ``Armpl_64ilp_mp``    (armpl, threaded code, ilp64 model)
  * ``Armpl_64ilp``       (armpl, sequential code, ilp64 model)
  * ``ACML``              (TODO)
  * ``ACML_MP``           (TODO)
  * ``ACML_GPU`           (TODO)

``THREADS``
  If set, checks only the thread version of libraries.

Result Variables
^^^^^^^^^^^^^^^^

This module defines the following variables:

``VENDOR``
  It contains the name of the vendor found if not specified
``VENDOR_FOUND``
  library implementing the vendor interface is found
``VENDOR_LIBRARIES``
  uncached list of vendor libraries (using full path name)
``BLA_VENDOR``
  Set the blas vendor
``BLAS_FOUND``
  library implementing the blas interface is found
``BLAS_LIBRARIES``
  uncached list of blas libraries (using full path name)
``LAPACK_FOUND``
  library implementing the lapack interface is found
``LAPACK_LIBRARIES``
  uncached list of lapack libraries (using full path name)
``FFTW_FOUND``
  library implementing the fftw interface is found
``FFTW_LIBRARIES``
  uncached list of fftw libraries (using full path name)
``FFTW_INCLUDE_DIRS``
  uncached list of include directory of fftw (using full path name)

.. note::

  C, CXX or Fortran must be enabled to detect a vendor library.

  For example, to use Intel MKL libraries and/or Intel compiler:

  .. code-block:: cmake

    set(VENDOR Intel10_64lp)
    find_package(VENDOR)

Hints
^^^^^

Set the ``MKLROOT`` environment variable to a directory that contains an MKL
installation, or add the directory to the dynamic library loader environment
variable for your platform (``LIB``, ``DYLD_LIBRARY_PATH`` or
``LD_LIBRARY_PATH``).

#]=======================================================================]

# Check the language being used
if(NOT (CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED OR CMAKE_Fortran_COMPILER_LOADED))
  if(VENDOR_FIND_REQUIRED)
    message(FATAL_ERROR "FindVENDOR requires Fortran, C, or C++ to be enabled.")
  else()
    message(STATUS "Looking for VENDOR... - NOT found (Unsupported languages)")
    return()
  endif()
endif()

if(NOT VENDOR_FIND_QUIETLY)
  message(STATUS "Looking for vendor")
endif()

if(CMAKE_Fortran_COMPILER_LOADED)
  include(CheckFortranFunctionExists)
else()
  include(CheckFunctionExists)
endif()
include(CMakePushCheckState)
include(FindPackageHandleStandardArgs)
cmake_push_check_state()
set(CMAKE_REQUIRED_QUIET ${VENDOR_FIND_QUIETLY})

set(_vendor_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})
if(VENDOR_STATIC)
  if(WIN32)
    set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
  else()
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
  endif()
else()
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # for ubuntu's libblas3gf and liblapack3gf packages
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
  endif()
endif()

set(VENDOR_LOOK_FOR_THREADS OFF)
if(VENDOR_FIND_COMPONENTS)
  foreach(component ${VENDOR_FIND_COMPONENTS})
    if (${component} STREQUAL "THREADS")
      # means we look for the Threads version of VENDOR
      set(VENDOR_LOOK_FOR_THREADS ON)
    endif()
  endforeach()
endif()

macro(CHECK_VENDOR_LIBRARIES LIBRARIES _prefix _name _flags _list _threadlibs _addlibdir _subdirs)
  # This macro checks for the existence of the combination of fortran libraries
  # given by _list.  If the combination is found, this macro checks (using the
  # Check_Fortran_Function_Exists macro) whether can link against that library
  # combination using the name of a routine given by _name using the linker
  # flags given by _flags.  If the combination of libraries is found and passes
  # the link test, LIBRARIES is set to the list of complete library paths that
  # have been found.  Otherwise, LIBRARIES is set to FALSE.

  # N.B. _prefix is the prefix applied to the names of all cached variables that
  # are generated internally and marked advanced by this macro.
  # _addlibdir is a list of additional search paths. _subdirs is a list of path
  # suffixes to be used by find_library().

  set(_libraries_work TRUE)
  set(${LIBRARIES})
  set(_combined_name)

  set(_extaddlibdir)
  if(WIN32)
    list(APPEND _extaddlibdir $ENV{LIB})
  elseif(APPLE)
    list(APPEND _extaddlibdir $ENV{DYLD_LIBRARY_PATH})
  else()
    list(APPEND _extaddlibdir $ENV{LD_LIBRARY_PATH})
  endif()
  string(REPLACE ":" ";" _extaddlibdir ${_extaddlibdir})
  list(APPEND _extaddlibdir "${_addlibdir}")
  list(APPEND _extaddlibdir "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")

  foreach(_library ${_list})
    if(_library MATCHES "^-Wl,--(start|end)-group$")
      # Respect linker flags like --start/end-group (required by MKL)
      set(${LIBRARIES} ${${LIBRARIES}} "${_library}")
    else()
      set(_combined_name ${_combined_name}_${_library})
      if(NOT "${_threadlibs}" STREQUAL "")
        set(_combined_name ${_combined_name}_threadlibs)
      endif()
      if(_libraries_work)
        find_library(${_prefix}_${_library}_LIBRARY
          NAMES ${_library}
          PATHS ${_extaddlibdir}
          PATH_SUFFIXES ${_subdirs}
        )
        #message("DEBUG: find_library(${_library}) got ${${_prefix}_${_library}_LIBRARY}")
        mark_as_advanced(${_prefix}_${_library}_LIBRARY)
        set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
        set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
      endif()
    endif()
  endforeach()

  if(_libraries_work)
    # Test this combination of libraries.
    set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_threadlibs})
    if(CMAKE_Fortran_COMPILER_LOADED)
      check_fortran_function_exists("${_name}" ${_prefix}${_combined_name}_WORKS)
    else()
      check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
    endif()
    set(CMAKE_REQUIRED_LIBRARIES)
    set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  endif()

  if(_libraries_work)
    if("${_list}" STREQUAL "")
      set(${LIBRARIES} "${LIBRARIES}-PLACEHOLDER-FOR-EMPTY-LIBRARIES")
    else()
      set(${LIBRARIES} ${${LIBRARIES}} ${_threadlibs})
    endif()
  else()
    set(${LIBRARIES} FALSE)
  endif()
  #message("DEBUG: ${LIBRARIES} = ${${LIBRARIES}}")
endmacro()

macro(CHECK_INCLUDE_DIRS INCLUDE_DIRS _prefix _list _addincdir _subdirs)
  # This macro checks for the existence of the combination of include libraries
  # that containthe _list files, INCLUDE_DIRS is set to the list of complete include paths 
  # that have been found. Otherwise, INCLUDE_DIRS is set to FALSE.

  # N.B. _prefix is the prefix applied to the names of all cached variables that
  # are generated internally and marked advanced by this macro.
  # _addincdir is a list of additional search paths. _subdirs is a list of path
  # suffixes to be used by find_path().
  
  set(${INCLUDE_DIRS})

  list(APPEND _extaincdir $ENV{CPATH})
  list(APPEND _extaincdir $ENV{C_INCLUDE_PATH})
  list(APPEND _extaincdir $ENV{CPLUS_INCLUDE_PATH})

  string(REPLACE ":" ";" _extaincdir ${_extaincdir})

  find_path(${_prefix}_INCDIR
    NAMES ${_list}
    PATHS ${_extaincdir} ${_addincdir}
    PATH_SUFFIXES ${_subdirs})
  #message("DEBUG: find_path( ${${_prefix}_INCDIR} ${_list} ${_extaincdir} ${_subdirs} )")

  if(${_prefix}_INCDIR)
    mark_as_advanced(${_prefix}_INCDIR)
    set(${INCLUDE_DIRS} ${${_prefix}_INCDIR})
  else()
    set(${INCLUDE_DIRS} FALSE)
  endif()
  #message("DEBUG: ${INCLUDE_DIRS} = ${${INCLUDE_DIRS}}")
endmacro()

set(VENDOR_LIBRARIES)
set(BLAS_LIBRARIES)
set(LAPACK_LIBRARIES)
set(FFTW_LIBRARIES)
set(FFTW_INCLUDE_DIRS)

if(NOT $ENV{VENDOR} STREQUAL "")
  set(VENDOR $ENV{VENDOR})
else()
  if(NOT VENDOR)
    set(VENDOR "All")
  endif()
endif()

# Intel MKL 10+ library?
if(VENDOR MATCHES "Intel" OR VENDOR STREQUAL "All")
  # Force VENDOR_LOOK_FOR_THREADS to be OFF if Intel VENDOR requires sequential code
  if(VENDOR MATCHES "_seq")
    set(VENDOR_LOOK_FOR_THREADS OFF)
  elseif(VENDOR MATCHES "Intel")
    set(VENDOR_LOOK_FOR_THREADS ON)
  endif()
  # System-specific settings
  if(WIN32)
    if(VENDOR_STATIC)
      set(VENDOR_mkl_DLL_SUFFIX "")
    else()
      set(VENDOR_mkl_DLL_SUFFIX "_dll")
    endif()
  else()
    if(VENDOR_STATIC)
      set(VENDOR_mkl_START_GROUP "-Wl,--start-group")
      set(VENDOR_mkl_END_GROUP "-Wl,--end-group")
    else()
      set(VENDOR_mkl_START_GROUP "")
      set(VENDOR_mkl_END_GROUP "")
    endif()
    # Switch to GNU Fortran support layer if needed (but not on Apple, where MKL does not provide it)
    if(CMAKE_Fortran_COMPILER_LOADED AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND NOT APPLE)
        set(VENDOR_mkl_INTFACE "gf")
        set(VENDOR_mkl_THREADING "gnu")
        set(VENDOR_mkl_OMP "gomp")
    else()
        set(VENDOR_mkl_INTFACE "intel")
        set(VENDOR_mkl_THREADING "intel")
        set(VENDOR_mkl_OMP "iomp5")
    endif()
    set(VENDOR_mkl_LM "-lm")
    set(VENDOR_mkl_LDL "-ldl")
  endif()

  if(VENDOR_FIND_QUIETLY OR NOT VENDOR_FIND_REQUIRED)
    find_package(Threads)
  else()
    find_package(Threads REQUIRED)
  endif()

  if(VENDOR MATCHES "_64ilp")
    set(VENDOR_mkl_ILP_MODE "ilp64")
  else()
    set(VENDOR_mkl_ILP_MODE "lp64")
  endif()

  set(VENDOR_SEARCH_LIBS "")

  set(VENDOR_mkl_SEARCH_SYMBOL sgemm)
  set(_LIBRARIES VENDOR_LIBRARIES)
  if(WIN32)
    # Find the main file (32-bit or 64-bit)
    set(VENDOR_SEARCH_LIBS_WIN_MAIN "")
    if(VENDOR STREQUAL "Intel10_32" OR VENDOR STREQUAL "All")
      list(APPEND VENDOR_SEARCH_LIBS_WIN_MAIN
        "mkl_intel_c${VENDOR_mkl_DLL_SUFFIX}")
    endif()
    if(VENDOR MATCHES "^Intel10_64i?lp" OR VENDOR STREQUAL "All")
      list(APPEND VENDOR_SEARCH_LIBS_WIN_MAIN
        "mkl_intel_${VENDOR_mkl_ILP_MODE}${VENDOR_mkl_DLL_SUFFIX}")
    endif()

    # Add threading/sequential libs
    set(VENDOR_SEARCH_LIBS_WIN_THREAD "")
    if(VENDOR_LOOK_FOR_THREADS)
      if(VENDOR MATCHES "^Intel10_64i?lp$" OR VENDOR STREQUAL "All")
        # old version
        list(APPEND VENDOR_SEARCH_LIBS_WIN_THREAD
          "libguide40 mkl_intel_thread${VENDOR_mkl_DLL_SUFFIX}")
        # mkl >= 10.3
        list(APPEND VENDOR_SEARCH_LIBS_WIN_THREAD
          "libiomp5md mkl_intel_thread${VENDOR_mkl_DLL_SUFFIX}")
      endif()
    else()
      if(VENDOR MATCHES "^Intel10_64i?lp_seq$" OR VENDOR STREQUAL "All")
        list(APPEND VENDOR_SEARCH_LIBS_WIN_THREAD
          "mkl_sequential${VENDOR_mkl_DLL_SUFFIX}")
      endif()
    endif()

    # Cartesian product of the above
    foreach(MAIN ${VENDOR_SEARCH_LIBS_WIN_MAIN})
      foreach(THREAD ${VENDOR_SEARCH_LIBS_WIN_THREAD})
        list(APPEND VENDOR_SEARCH_LIBS
          "${MAIN} ${THREAD} mkl_core${VENDOR_mkl_DLL_SUFFIX}")
      endforeach()
    endforeach()
  else()
    if(VENDOR STREQUAL "Intel10_32" OR VENDOR STREQUAL "All")
      set(SEARCH_LIBS "mkl_${VENDOR_mkl_INTFACE}")
      if(VENDOR_LOOK_FOR_THREADS)
        string(CONCAT SEARCH_LIBS " mkl_${VENDOR_mkl_THREADING}_thread")
      endif()
      string(CONCAT SEARCH_LIBS " mkl_core guide")
      list(APPEND VENDOR_SEARCH_LIBS ${SEARCH_LIBS})

      # mkl >= 10.3
      set(SEARCH_LIBS "${VENDOR_mkl_START_GROUP} mkl_${VENDOR_mkl_INTFACE}")
      if(VENDOR_LOOK_FOR_THREADS)
        string(CONCAT SEARCH_LIBS " mkl_${VENDOR_mkl_THREADING}_thread")
      endif()
      string(CONCAT SEARCH_LIBS " mkl_core ${VENDOR_mkl_END_GROUP}")
      if(VENDOR_LOOK_FOR_THREADS)
        string(CONCAT SEARCH_LIBS " ${VENDOR_mkl_OMP}")
      endif()
      list(APPEND VENDOR_SEARCH_LIBS ${SEARCH_LIBS})
    endif()
    if(VENDOR_LOOK_FOR_THREADS)
      if(VENDOR MATCHES "^Intel10_64i?lp$" OR VENDOR STREQUAL "All")
        # old version
        list(APPEND VENDOR_SEARCH_LIBS
          "mkl_${VENDOR_mkl_INTFACE}_${VENDOR_mkl_ILP_MODE} mkl_${VENDOR_mkl_THREADING}_thread mkl_core guide")

        # mkl >= 10.3
        list(APPEND VENDOR_SEARCH_LIBS
          "${VENDOR_mkl_START_GROUP} mkl_${VENDOR_mkl_INTFACE}_${VENDOR_mkl_ILP_MODE} mkl_${VENDOR_mkl_THREADING}_thread mkl_core ${VENDOR_mkl_END_GROUP} ${VENDOR_mkl_OMP}")
      endif()
    else()
      if(VENDOR MATCHES "^Intel10_64i?lp_seq$" OR VENDOR STREQUAL "All")
        list(APPEND VENDOR_SEARCH_LIBS
          "${VENDOR_mkl_START_GROUP} mkl_${VENDOR_mkl_INTFACE}_${VENDOR_mkl_ILP_MODE} mkl_sequential mkl_core ${VENDOR_mkl_END_GROUP}")
      endif()
    endif()

    #older vesions of intel mkl libs
    if(VENDOR STREQUAL "Intel" OR VENDOR STREQUAL "All")
      list(APPEND VENDOR_SEARCH_LIBS
        "mkl")
      list(APPEND VENDOR_SEARCH_LIBS
        "mkl_ia32")
      list(APPEND VENDOR_SEARCH_LIBS
        "mkl_em64t")
    endif()
  endif()

  if(VENDOR MATCHES "^Intel10_64_dyn$" OR VENDOR STREQUAL "All")
    # mkl >= 10.3 with single dynamic library
    list(APPEND VENDOR_SEARCH_LIBS
      "mkl_rt")
  endif()

  # MKL uses a multitude of partially platform-specific subdirectories:
  if(VENDOR STREQUAL "Intel10_32")
    set(VENDOR_mkl_ARCH_NAME "ia32")
  else()
    set(VENDOR_mkl_ARCH_NAME "intel64")
  endif()
  if(WIN32)
    set(VENDOR_mkl_OS_NAME "win")
  elseif(APPLE)
    set(VENDOR_mkl_OS_NAME "mac")
  else()
    set(VENDOR_mkl_OS_NAME "lin")
  endif()
  if(DEFINED ENV{MKLROOT})
    set(VENDOR_mkl_MKLROOT "$ENV{MKLROOT}")
    # If MKLROOT points to the subdirectory 'mkl', use the parent directory instead
    # so we can better detect other relevant libraries in 'compiler' or 'tbb':
    get_filename_component(VENDOR_mkl_MKLROOT_LAST_DIR "${VENDOR_mkl_MKLROOT}" NAME)
    if(VENDOR_mkl_MKLROOT_LAST_DIR STREQUAL "mkl")
        get_filename_component(VENDOR_mkl_MKLROOT "${VENDOR_mkl_MKLROOT}" DIRECTORY)
    endif()
  endif()
  set(VENDOR_mkl_LIB_PATH_SUFFIXES
      "compiler/lib" "compiler/lib/${VENDOR_mkl_ARCH_NAME}_${VENDOR_mkl_OS_NAME}"
      "mkl/lib" "mkl/lib/${VENDOR_mkl_ARCH_NAME}_${VENDOR_mkl_OS_NAME}"
      "lib/${VENDOR_mkl_ARCH_NAME}_${VENDOR_mkl_OS_NAME}")

  foreach(IT ${VENDOR_SEARCH_LIBS})
    string(REPLACE " " ";" SEARCH_LIBS ${IT})
    if(NOT ${_LIBRARIES})
      check_vendor_libraries(
        ${_LIBRARIES}
        VENDOR
        ${VENDOR_mkl_SEARCH_SYMBOL}
        ""
        "${SEARCH_LIBS}"
        "${CMAKE_THREAD_LIBS_INIT};${VENDOR_mkl_LM};${VENDOR_mkl_LDL}"
        "${VENDOR_mkl_MKLROOT}"
        "${VENDOR_mkl_LIB_PATH_SUFFIXES}"
        )
    endif()
  endforeach()
  set(BLAS_LIBRARIES ${VENDOR_LIBRARIES})
  set(LAPACK_LIBRARIES ${VENDOR_LIBRARIES})
  set(FFTW_LIBRARIES ${VENDOR_LIBRARIES})

  if(VENDOR MATCHES "All" AND VENDOR_LIBRARIES)
    if(VENDOR_LIBRARIES MATCHES "ilp64")
      if(VENDOR_LIBRARIES MATCHES "sequential")
        set(VENDOR "Intel10_64ilp_seq")
      else()
        set(VENDOR "Intel10_64ilp")
      endif()
    else()
      if(VENDOR_LIBRARIES MATCHES "sequential")
        set(VENDOR "Intel10_64lp_seq")
      else()
        set(VENDOR "Intel10_64lp")
      endif()
    endif()
  endif()

  unset(VENDOR_mkl_ILP_MODE)
  unset(VENDOR_mkl_INTFACE)
  unset(VENDOR_mkl_THREADING)
  unset(VENDOR_mkl_OMP)
  unset(VENDOR_mkl_DLL_SUFFIX)
  unset(VENDOR_mkl_LM)
  unset(VENDOR_mkl_LDL)
  unset(VENDOR_mkl_MKLROOT)
  unset(VENDOR_mkl_MKLROOT_LAST_DIR)
  unset(VENDOR_mkl_ARCH_NAME)
  unset(VENDOR_mkl_OS_NAME)
  unset(VENDOR_mkl_LIB_PATH_SUFFIXES)
endif()

# ARMPL library?
if(VENDOR MATCHES "Armpl" OR VENDOR STREQUAL "All")
  if(CMAKE_C_COMPILER_LOADED OR CMAKE_CXX_COMPILER_LOADED)
    # System-specific settings
    set(VENDOR_armpl_LAMATH "-lamath")
    set(VENDOR_armpl_LM "-lm")

    set(VENDOR_armpl_SEARCH)
    if(VENDOR MATCHES "_64ilp")
      set(VENDOR_armpl_SEARCH "armpl_ilp64")
    else()
      set(VENDOR_armpl_SEARCH "armpl_lp64")
    endif()

    # Force VENDOR_LOOK_FOR_THREADS to be ON if ARMPL VENDOR requires Threads
    if(VENDOR MATCHES "_mp")
      set(VENDOR_LOOK_FOR_THREADS ON)
    elseif(VENDOR MATCHES "Armpl")
      set(VENDOR_LOOK_FOR_THREADS OFF)
    endif()

    if(VENDOR_LOOK_FOR_THREADS)
      if(NOT Threads_FOUND)
        if(VENDOR_FIND_QUIETLY OR NOT VENDOR_FIND_REQUIRED)
          find_package(Threads)
        else()
          find_package(Threads REQUIRED)
        endif()
      endif()
      string(CONCAT VENDOR_armpl_SEARCH ${VENDOR_armpl_SEARCH} "_mp")
      string(CONCAT VENDOR_armpl_SEARCH ${VENDOR_armpl_SEARCH} ";omp")
    endif()

    if(NOT VENDOR_LIBRARIES)
      check_vendor_libraries(
        VENDOR_LIBRARIES
        VENDOR
        sgemm
        ""
        "${VENDOR_armpl_SEARCH}"
        "${VENDOR_armpl_LAMATH};${VENDOR_armpl_LM}"
        ""
        ""
        )
        set(BLAS_LIBRARIES ${VENDOR_LIBRARIES})
        set(LAPACK_LIBRARIES ${VENDOR_LIBRARIES})
        set(FFTW_LIBRARIES ${VENDOR_LIBRARIES})
    endif()

    if(VENDOR MATCHES "All" AND VENDOR_LIBRARIES)
      if(VENDOR_LIBRARIES MATCHES "ilp64")
        if(VENDOR_LIBRARIES MATCHES "_mp")
          set(VENDOR "Armpl_64ilp_mp")
        else()
          set(VENDOR "Armpl_64ilp")
        endif()
      else()
        if(VENDOR_LIBRARIES MATCHES "_mp")
          set(VENDOR "Armpl_64lp_mp")
        else()
          set(VENDOR "Armpl_64lp")
        endif()
      endif()
    endif()

    unset(VENDOR_armpl_SEARCH)
    unset(VENDOR_armpl_LAMATH)
    unset(VENDOR_armpl_LM)
  endif()
endif()

## ACML library?
#if(VENDOR MATCHES "ACML" OR VENDOR STREQUAL "All")
#  if(((VENDOR STREQUAL "ACML") AND (NOT VENDOR_ACML_LIB_DIRS)) OR
#    ((VENDOR STREQUAL "ACML_MP") AND (NOT VENDOR_ACML_MP_LIB_DIRS)) OR
#    ((VENDOR STREQUAL "ACML_GPU") AND (NOT VENDOR_ACML_GPU_LIB_DIRS)))
#    # try to find acml in "standard" paths
#    if(WIN32)
#      file(GLOB _ACML_ROOT "C:/AMD/acml*/ACML-EULA.txt")
#    else()
#      file(GLOB _ACML_ROOT "/opt/acml*/ACML-EULA.txt")
#    endif()
#    if(WIN32)
#      file(GLOB _ACML_GPU_ROOT "C:/AMD/acml*/GPGPUexamples")
#    else()
#      file(GLOB _ACML_GPU_ROOT "/opt/acml*/GPGPUexamples")
#    endif()
#    list(GET _ACML_ROOT 0 _ACML_ROOT)
#    list(GET _ACML_GPU_ROOT 0 _ACML_GPU_ROOT)
#    if(_ACML_ROOT)
#      get_filename_component(_ACML_ROOT ${_ACML_ROOT} PATH)
#      if(SIZEOF_INTEGER EQUAL 8)
#        set(_ACML_PATH_SUFFIX "_int64")
#      else()
#        set(_ACML_PATH_SUFFIX "")
#      endif()
#      if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
#        set(_ACML_COMPILER32 "ifort32")
#        set(_ACML_COMPILER64 "ifort64")
#      elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro")
#        set(_ACML_COMPILER32 "sun32")
#        set(_ACML_COMPILER64 "sun64")
#      elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
#        set(_ACML_COMPILER32 "pgi32")
#        if(WIN32)
#          set(_ACML_COMPILER64 "win64")
#        else()
#          set(_ACML_COMPILER64 "pgi64")
#        endif()
#      elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Open64")
#        # 32 bit builds not supported on Open64 but for code simplicity
#        # We'll just use the same directory twice
#        set(_ACML_COMPILER32 "open64_64")
#        set(_ACML_COMPILER64 "open64_64")
#      elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
#        set(_ACML_COMPILER32 "nag32")
#        set(_ACML_COMPILER64 "nag64")
#      else()
#        set(_ACML_COMPILER32 "gfortran32")
#        set(_ACML_COMPILER64 "gfortran64")
#      endif()
#
#      if(VENDOR STREQUAL "ACML_MP")
#        set(_ACML_MP_LIB_DIRS
#          "${_ACML_ROOT}/${_ACML_COMPILER32}_mp${_ACML_PATH_SUFFIX}/lib"
#          "${_ACML_ROOT}/${_ACML_COMPILER64}_mp${_ACML_PATH_SUFFIX}/lib")
#      else()
#        set(_ACML_LIB_DIRS
#          "${_ACML_ROOT}/${_ACML_COMPILER32}${_ACML_PATH_SUFFIX}/lib"
#          "${_ACML_ROOT}/${_ACML_COMPILER64}${_ACML_PATH_SUFFIX}/lib")
#      endif()
#    endif()
#  elseif(VENDOR_${VENDOR}_LIB_DIRS)
#    set(_${VENDOR}_LIB_DIRS ${VENDOR_${VENDOR}_LIB_DIRS})
#  endif()
#
#  if(VENDOR STREQUAL "ACML_MP" OR VENDOR_LOOK_FOR_THREADS)
#    foreach(VENDOR_ACML_MP_LIB_DIRS ${_ACML_MP_LIB_DIRS})
#      check_vendor_libraries(
#        VENDOR_LIBRARIES
#        VENDOR
#        sgemm
#        "" "acml_mp;acml_mv" "" ${VENDOR_ACML_MP_LIB_DIRS} ""
#        )
#      set(BLAS_LIBRARIES ${VENDOR_LIBRARIES})
#      set(LAPACK_LIBRARIES ${VENDOR_LIBRARIES})
#      set(FFTW_LIBRARIES ${VENDOR_LIBRARIES})
#      if(VENDOR_LIBRARIES)
#        break()
#      endif()
#    endforeach()
#    if(VENDOR MATCHES "All" AND VENDOR_LIBRARIES)
#      set(VENDOR "ACML_MP")
#    endif()
#  elseif(VENDOR STREQUAL "ACML_GPU")
#    foreach(VENDOR_ACML_GPU_LIB_DIRS ${_ACML_GPU_LIB_DIRS})
#      check_vendor_libraries(
#        VENDOR_LIBRARIES
#        VENDOR
#        sgemm
#        "" "acml;acml_mv;CALBLAS" "" ${VENDOR_ACML_GPU_LIB_DIRS} ""
#        )
#      set(BLAS_LIBRARIES ${VENDOR_LIBRARIES})
#      set(LAPACK_LIBRARIES ${VENDOR_LIBRARIES})
#      set(FFTW_LIBRARIES ${VENDOR_LIBRARIES})
#      if(VENDOR_LIBRARIES)
#        break()
#      endif()
#    endforeach()
#    if(VENDOR MATCHES "All" AND VENDOR_LIBRARIES)
#      set(VENDOR "ACML_GPU")
#    endif()
#  else()
#    foreach(VENDOR_ACML_LIB_DIRS ${_ACML_LIB_DIRS})
#      check_vendor_libraries(
#        VENDOR_LIBRARIES
#        VENDOR
#        sgemm
#        "" "acml;acml_mv" "" ${VENDOR_ACML_LIB_DIRS} ""
#        )
#      set(BLAS_LIBRARIES ${VENDOR_LIBRARIES})
#      set(LAPACK_LIBRARIES ${VENDOR_LIBRARIES})
#      set(FFTW_LIBRARIES ${VENDOR_LIBRARIES})
#      if(VENDOR_LIBRARIES)
#        break()
#      endif()
#    endforeach()
#    if(VENDOR MATCHES "All" AND VENDOR_LIBRARIES)
#      set(VENDOR "ACML")
#    endif()
#  endif()
#
#  # Either acml or acml_mp should be in LD_LIBRARY_PATH but not both
#  if(NOT VENDOR_LIBRARIES)
#    check_vendor_libraries(
#      VENDOR_LIBRARIES
#      VENDOR
#      sgemm
#      ""
#      "acml;acml_mv"
#      ""
#      ""
#      ""
#      )
#    set(BLAS_LIBRARIES ${VENDOR_LIBRARIES})
#    set(LAPACK_LIBRARIES ${VENDOR_LIBRARIES})
#    set(FFTW_LIBRARIES ${VENDOR_LIBRARIES})
#    if(VENDOR MATCHES "All" AND VENDOR_LIBRARIES)
#      set(VENDOR "ACML")
#    endif()
#  endif()
#  if(NOT VENDOR_LIBRARIES)
#    check_vendor_libraries(
#      VENDOR_LIBRARIES
#      VENDOR
#      sgemm
#      ""
#      "acml_mp;acml_mv"
#      ""
#      ""
#      ""
#      )
#    set(BLAS_LIBRARIES ${VENDOR_LIBRARIES})
#    set(LAPACK_LIBRARIES ${VENDOR_LIBRARIES})
#    set(FFTW_LIBRARIES ${VENDOR_LIBRARIES})
#    if(VENDOR MATCHES "All" AND VENDOR_LIBRARIES)
#      set(VENDOR "ACML_MP")
#    endif()
#  endif()
#  if(NOT VENDOR_LIBRARIES)
#    check_vendor_libraries(
#      VENDOR_LIBRARIES
#      VENDOR
#      sgemm
#      ""
#      "acml;acml_mv;CALBLAS"
#      ""
#      ""
#      ""
#      )
#    set(BLAS_LIBRARIES ${VENDOR_LIBRARIES})
#    set(LAPACK_LIBRARIES ${VENDOR_LIBRARIES})
#    set(FFTW_LIBRARIES ${VENDOR_LIBRARIES})
#    if(VENDOR MATCHES "All" AND VENDOR_LIBRARIES)
#      set(VENDOR "ACML_GPU")
#    endif()
#  endif()
#endif() # ACML

if(NOT FFTW_INCLUDE_DIRS)
  set(_libdirs)

  foreach(_elem ${VENDOR_LIBRARIES})
    if(_elem MATCHES "^\/")
      get_filename_component(_libdir ${_elem} PATH)
      get_filename_component(_libdir ${_libdir} PATH)
      list(APPEND _libdirs "${_libdir}/include")
      get_filename_component(_libdir ${_libdir} PATH)
      list(APPEND _libdirs "${_libdir}/include")
    endif()
  endforeach()

  if(VENDOR MATCHES "Intel")
    list(APPEND _libdirs 
      "$ENV{MKL_INC}"
      "$ENV{MKL_INCLUDE}"
      "$ENV{MKL_INCLUDES}"
      "$ENV{MKLROOT}/include"
      "$ENV{MKL_HOME}/include"
      "$ENV{MKL_ROOT}/include"
      "$ENV{MKL_PATH}/include"
      "$ENV{MKL_DIR}/include")
  elseif(VENDOR MATCHES "Armpl")
    list(APPEND _libdirs 
      "$ENV{ARMPL_INC}"
      "$ENV{ARMPL_INCLUDE}"
      "$ENV{ARMPL_INCLUDES}"
      "$ENV{ARMPL_HOME}/include"
      "$ENV{ARMPL_ROOT}/include"
      "$ENV{ARMPL_PATH}/include"
      "$ENV{ARMPL_DIR}/include")
  elseif(VENDOR MATCHES "ACML")
    # TODO
  endif()

  list(REMOVE_DUPLICATES _libdirs)

  check_include_dirs(
    FFTW_INCLUDE_DIRS
    FFTW
    "fftw3.h;fftw3.f"
    "${_libdirs}"
    "fftw")

  if(VENDOR MATCHES "Intel")
    check_include_dirs(
      FFTW_MKL_DFTI_F90_DIR
      MKL_DFTI
      "mkl_dfti.f90"
      "${_libdirs}"
      "")
    list(APPEND FFTW_INCLUDE_DIRS ${FFTW_MKL_DFTI_F90_DIR})
  endif()
endif()

# On compilers that implicitly link VENDOR (such as ftn, cc, and CC on Cray HPC machines)
# we used a placeholder for empty VENDOR_LIBRARIES to get through our logic above.
if(VENDOR_LIBRARIES STREQUAL "VENDOR_LIBRARIES-PLACEHOLDER-FOR-EMPTY-LIBRARIES")
  set(VENDOR_LIBRARIES "")
endif()

cmake_pop_check_state()
set(CMAKE_FIND_LIBRARY_SUFFIXES ${_vendor_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

if(VENDOR_LIBRARIES AND FFTW_INCLUDE_DIRS)
  set(VENDOR_FOUND TRUE)
  set(BLAS_FOUND TRUE)
  set(LAPACK_FOUND TRUE)
  set(FFTW_FOUND TRUE)
  set(BLA_VENDOR ${VENDOR})
else()
  set(VENDOR_FOUND FALSE)
  set(BLAS_FOUND FALSE)
  set(LAPACK_FOUND FALSE)
  set(FFTW_FOUND FALSE)

  if(VENDOR_FIND_REQUIRED)
    message(FATAL_ERROR "Vendor not found!")
  else()
    if(NOT VENDOR_FIND_QUIETLY)
      message(STATUS "Vendor not found!")
    endif()
  endif()
endif()