#[=======================================================================[.rst:
FindFoX
----------

Find FoX library

This module finds an installed FoX library.

Input Variables
^^^^^^^^^^^^^^^

The following variables may be set to influence this module's behavior:

``FOX_ROOT``
  The path to the installation folder of FoX library

Imported targets
^^^^^^^^^^^^^^^^

This module defines the following :prop_tgt:`IMPORTED` targets:

``FoX::FoX``
``FoX::DOM``
``FoX::SAX``
``FoX::WXML``
``FoX::Common``
``FoX::Utils``
``FoX::FSys``

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project:

``FOX_FOUND``
  FoX is found
``FOX_LIBRARIES``
  the libraries needed to use FoX.
``FOX_INCLUDE_DIRS``
  where to find modules and headers for FoX
``FOX_DOM_LIB``
  DOM interface of FoX
``FOX_SAX_LIB``
  SAX interface of FoX
``FOX_WXML_LIB``
  VoiceXML interface of FoX
``FOX_COMMON_LIB``
  Interface for common functions of FoX
``FOX_UTILS_LIB``
  Interface for util functions of FoX
``FOX_FSYS_LIB``
  Interface for file system functions of FoX

#]=======================================================================]

find_library(
    FOX_DOM_LIB
    NAMES "FoX_dom"
    HINTS ${FOX_ROOT}
    PATH_SUFFIXES "lib")

find_library(
    FOX_SAX_LIB
    NAMES "FoX_sax"
    HINTS ${FOX_ROOT}
    PATH_SUFFIXES "lib")

find_library(
    FOX_WXML_LIB
    NAMES "FoX_wxml"
    HINTS ${FOX_ROOT}
    PATH_SUFFIXES "lib")

find_library(
    FOX_COMMON_LIB
    NAMES "FoX_common"
    HINTS ${FOX_ROOT}
    PATH_SUFFIXES "lib")

find_library(
    FOX_UTILS_LIB
    NAMES "FoX_utils"
    HINTS ${FOX_ROOT}
    PATH_SUFFIXES "lib")

find_library(
    FOX_FSYS_LIB
    NAMES "FoX_fsys"
    HINTS ${FOX_ROOT}
    PATH_SUFFIXES "lib")

set(FOX_LIBRARIES 
    ${FOX_DOM_LIB}
    ${FOX_SAX_LIB}
    ${FOX_WXML_LIB}
    ${FOX_COMMON_LIB}
    ${FOX_UTILS_LIB}
    ${FOX_FSYS_LIB})

find_path(
    FOX_INCLUDE_DIRS
    NAMES "m_common_io.mod"
    HINTS ${FOX_ROOT}
    PATH_SUFFIXES
        "include"
        "finclude")

find_package_handle_standard_args(FOX
  REQUIRED_VARS
    FOX_LIBRARIES
    FOX_DOM_LIB
    FOX_SAX_LIB
    FOX_WXML_LIB
    FOX_COMMON_LIB
    FOX_UTILS_LIB
    FOX_FSYS_LIB
    FOX_INCLUDE_DIRS)

if(FOX_FOUND)
  add_library(FoX::FoX INTERFACE IMPORTED)
  add_library(FoX::DOM INTERFACE IMPORTED)
  add_library(FoX::SAX INTERFACE IMPORTED)
  add_library(FoX::WXML INTERFACE IMPORTED)
  add_library(FoX::Common INTERFACE IMPORTED)
  add_library(FoX::Utils INTERFACE IMPORTED)
  add_library(FoX::FSys INTERFACE IMPORTED)

  target_link_libraries(FoX::FoX
    INTERFACE 
      FoX::DOM
      FoX::SAX
      FoX::WXML
      FoX::Common
      FoX::Utils
      FoX::FSys)

  target_link_libraries(FoX::DOM INTERFACE ${FOX_DOM_LIB})
  target_link_libraries(FoX::SAX INTERFACE ${FOX_SAX_LIB})
  target_link_libraries(FoX::WXML INTERFACE ${FOX_WXML_LIB})
  target_link_libraries(FoX::Common INTERFACE ${FOX_COMMON_LIB})
  target_link_libraries(FoX::Utils INTERFACE ${FOX_UTILS_LIB})
  target_link_libraries(FoX::FSys INTERFACE ${FOX_FSYS_LIB})

  target_include_directories(FoX::DOM INTERFACE ${FOX_INCLUDE_DIRS})
  target_include_directories(FoX::SAX INTERFACE ${FOX_INCLUDE_DIRS})
  target_include_directories(FoX::WXML INTERFACE ${FOX_INCLUDE_DIRS})
  target_include_directories(FoX::Common INTERFACE ${FOX_INCLUDE_DIRS})
  target_include_directories(FoX::Utils INTERFACE ${FOX_INCLUDE_DIRS})
  target_include_directories(FoX::FSys INTERFACE ${FOX_INCLUDE_DIRS})
endif()
