#[=======================================================================[.rst:
FindWannier90
----------

Find Wannier90 library

This module finds an installed Wannier90 library.

Input Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``WANNIER90_ROOT``
  The path to the installation folder of Wannier90 library

Imported targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` targets:

``Wannier90::Wannier90``

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project:

``WANNIER90_FOUND``
  Wannier90 is found
``WANNIER90_LIBRARIES``
  the libraries needed to use Wannier90.
``WANNIER90_INCLUDE_DIRS``
  where to find modules and headers for Wannier90

#]=======================================================================]

find_library(
  WANNIER90_LIBRARIES
  NAMES "wannier"
  HINTS ${WANNIER90_ROOT}
  PATH_SUFFIXES "lib")

find_path(
  WANNIER90_INCLUDE_DIRS
  NAMES "w90_io.mod"
  HINTS ${WANNIER90_ROOT}
  PATH_SUFFIXES 
      "include"
      "modules")

find_package_handle_standard_args(WANNIER90
  REQUIRED_VARS
    WANNIER90_LIBRARIES
    WANNIER90_INCLUDE_DIRS)

if(WANNIER90_FOUND)
  add_library(Wannier90::Wannier90 INTERFACE IMPORTED)
  target_link_libraries(Wannier90::Wannier90 
    INTERFACE ${WANNIER90_LIBRARIES})
  target_include_directories(Wannier90::Wannier90
    INTERFACE ${WANNIER90_INCLUDE_DIRS})
endif()
