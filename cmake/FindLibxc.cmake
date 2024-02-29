include(FindPackageHandleStandardArgs)
find_package(PkgConfig REQUIRED)

pkg_search_module(_LIBXC libxc>=${Libxc_FIND_VERSION})

find_library(LIBXC_LIBRARIES NAMES libxc.a
  PATH_SUFFIXES lib
  HINTS
  ${LIBXC_ROOT}
  ENV EBROOTLIBXC
  ENV LIBXCROOT
  ${_LIBXC_LIBRARY_DIRS}
  DOC "libxc libraries list")

find_library(LIBXC_LIBRARIES_F03 NAMES libxcf03.a
  PATH_SUFFIXES lib
  HINTS
  ${LIBXC_ROOT}
  ENV EBROOTLIBXC
  ENV LIBXCROOT
  ${_LIBXC_LIBRARY_DIRS}
  DOC "libxc libraries list")

find_path(LIBXC_INCLUDE_DIR NAMES xc.h
  PATH_SUFFIXES inc include
  HINTS
  ${LIBXC_ROOT}
  ${_LIBXC_INCLUDE_DIRS}
  ENV EBROOTLIBXC
  ENV LIBXCROOT)

find_path(LIBXC_INCLUDE_DIR_F03 NAMES xc_f03_lib_m.mod
  PATH_SUFFIXES inc include
  HINTS
  ${LIBXC_ROOT}
  ${_LIBXC_INCLUDE_DIRS}
  ENV EBROOTLIBXC
  ENV LIBXCROOT)

find_package_handle_standard_args(Libxc DEFAULT_MSG LIBXC_LIBRARIES LIBXC_INCLUDE_DIR)

if (${Libxc_FOUND} AND LIBXC_LIBRARIES_F03 AND LIBXC_INCLUDE_DIR_F03)
  if(_LIBXC_VERSION)
    set(Libxc_VERSION ${_LIBXC_VERSION})
  else()
    set(Libxc_VERSION ${Libxc_FIND_VERSION})
  endif()
  set(Libxc_INCLUDE_DIR ${LIBXC_INCLUDE_DIR})
  add_library(Libxc::xcf03 INTERFACE IMPORTED)
  set_target_properties(Libxc::xcf03 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${LIBXC_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${LIBXC_LIBRARIES_F03};${LIBXC_LIBRARIES}")
else()
  unset(Libxc_FOUND)
endif()

