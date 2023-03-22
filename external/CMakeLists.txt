###########################################################
# FoX
###########################################################
include(fox.cmake)

###########################################################
# WANNIER90
###########################################################
include(wannier90.cmake)

###########################################################
# MBD
###########################################################
include(mbd.cmake)

###########################################################
# DeviceXlib
###########################################################
include(devxlib.cmake)

###########################################################
# Environ
###########################################################
include(environ.cmake)

###########################################################
# D3Q
###########################################################
if("d3q" IN_LIST QE_ENABLE_PLUGINS)
  message(STATUS "Installing D3Q via submodule")
  qe_git_submodule_update(external/d3q)
  add_subdirectory(d3q)
endif()

###########################################################
# PW2QMCPACK
###########################################################
if("pw2qmcpack" IN_LIST QE_ENABLE_PLUGINS OR QE_ENABLE_PW2QMCPACK)
  message(STATUS "Installing pw2qmcpack via submodule")
  qe_git_submodule_update(external/pw2qmcpack)
  add_subdirectory(pw2qmcpack)
endif()

###########################################################
# QE-GIPAW
###########################################################
if("gipaw" IN_LIST QE_ENABLE_PLUGINS)
  if(QE_ENABLE_FOX)
    message(STATUS "Installing qe-gipaw via submodule")
    qe_git_submodule_update(external/qe-gipaw)
    add_subdirectory(qe-gipaw)
  else()
    message(FATAL_ERROR "QE-GIPAW requires -DQE_ENABLE_FOX=ON .")
  endif()
endif()
