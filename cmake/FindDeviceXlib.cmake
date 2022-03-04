#[=======================================================================[.rst:
FindDEVICEXLIB
----------

Find DEVICEXLIB library

This module finds an installed DeviceXlib library.

Input Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``DEVICEXLIB_ROOT``
  The path to the installation folder of DeviceXlib library

Imported targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` targets:

``DeviceXlib::DeviceXlib``

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project:

``DEVICEXLIB_FOUND``
  DeviceXlib is found
``DEVICEXLIB_LIBRARIES``
  the libraries needed to use DeviceXlib.
``DEVICEXLIB_INCLUDE_DIRS``
  where to find modules and headers for DeviceXlib

#]=======================================================================]

find_library(
  DEVICEXLIB_LIBRARIES
  NAMES "devXlib"
  HINTS ${DEVICEXLIB_ROOT}
  PATH_SUFFIXES 
      "lib"
      "src")

find_path(
  DEVICEXLIB_INCLUDE_DIRS
  NAMES "device_fbuff_m.mod"
  HINTS ${DEVICEXLIB_ROOT}
  PATH_SUFFIXES 
      "include"
      "src")

find_package_handle_standard_args(DEVICEXLIB
  REQUIRED_VARS
    DEVICEXLIB_LIBRARIES
    DEVICEXLIB_INCLUDE_DIRS)

if(DEVICEXLIB_FOUND)
  add_library(DeviceXlib::DeviceXlib INTERFACE IMPORTED)
  target_link_libraries(DeviceXlib::DeviceXlib
    INTERFACE ${DEVICEXLIB_LIBRARIES})
  target_include_directories(DeviceXlib::DeviceXlib
    INTERFACE ${DEVICEXLIB_INCLUDE_DIRS})
endif()
