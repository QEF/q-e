# Check compiler version
if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 4.9)
  message(FATAL_ERROR "Requires GCC 4.9 or higher ")
endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
endif()
