#[=======================================================================[.rst:
FindMBD
----------

Find MBD library

This module finds an installed MBD library.

Input Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``MBD_ROOT``
  The path to the installation folder of MBD library

Imported targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` targets:

``MBD::MBD``

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project:

``MBD_FOUND``
  MBD is found
``MBD_LIBRARIES``
  the libraries needed to use MBD.
``MBD_INCLUDE_DIRS``
  where to find modules and headers for MBD

#]=======================================================================]

find_library(
    MBD_LIBRARIES
    NAMES "mbd"
    HINTS ${MBD_ROOT}
    PATH_SUFFIXES "lib")

find_path(
    MBD_INCLUDE_DIRS
    NAMES "mbd.mod"
    HINTS ${MBD_ROOT}
    PATH_SUFFIXES 
        "mbd"
        "include")

find_package_handle_standard_args(MBD
  REQUIRED_VARS
    MBD_LIBRARIES
    MBD_INCLUDE_DIRS)

if(MBD_FOUND)
  add_library(MBD::MBD INTERFACE IMPORTED)
  target_link_libraries(MBD::MBD 
    INTERFACE ${MBD_LIBRARIES})
  target_include_directories(MBD::MBD
    INTERFACE ${MBD_INCLUDE_DIRS})
endif()
