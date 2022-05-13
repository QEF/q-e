#[=======================================================================[.rst:
FindEnviron
----------

Find Environ library

This module finds an installed Environ library.

Input Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``ENVIRON_ROOT``
  The path to the installation folder of Environ library

Imported targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` targets:

``Environ::Environ``

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project:

``ENVIRON_FOUND``
  Environ is found
``ENVIRON_UTIL_LIB``
  Environ's version of QE's UtilXlib.
``ENVIRON_FFT_LIB``
  Environ's version of QE's FFTXlib.
``ENVIRON_SRC_LIB``
  Environ's source library.
``ENVIRON_INCLUDE_DIRS``
  where to find modules and headers for Environ

#]=======================================================================]

find_library(
  ENVIRON_UTIL_LIB
  NAMES "libenvutil.a"
  HINTS ${ENVIRON_ROOT}
  PATH_SUFFIXES
      "libs")

find_library(
  ENVIRON_FFT_LIB
  NAMES "libenvfft.a"
  HINTS ${ENVIRON_ROOT}
  PATH_SUFFIXES
  "libs")

find_library(
  ENVIRON_SRC_LIB
  NAMES "libenvsrc.a"
  HINTS ${ENVIRON_ROOT}
  PATH_SUFFIXES
  "libs")

find_path(
  ENVIRON_INCLUDE_SRC
  NAMES "environ_api.mod"
  HINTS ${ENVIRON_ROOT}
  PATH_SUFFIXES
      "src")

find_package_handle_standard_args(Environ
  REQUIRED_VARS
  ENVIRON_UTIL_LIB
  ENVIRON_FFT_LIB
  ENVIRON_SRC_LIB
  ENVIRON_INCLUDE_SRC)

if(ENVIRON_FOUND)
  add_library(Environ::Environ INTERFACE IMPORTED)
  target_link_libraries(Environ::Environ
    INTERFACE ${ENVIRON_SRC_LIB})
  target_link_libraries(Environ::Environ
    INTERFACE ${ENVIRON_FFT_LIB})
  target_link_libraries(Environ::Environ
    INTERFACE ${ENVIRON_UTIL_LIB})
  target_include_directories(Environ::Environ
    INTERFACE ${ENVIRON_INCLUDE_SRC})
endif()
