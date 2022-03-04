# Check compiler version
if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 4.9)
  message(FATAL_ERROR "Requires GCC 4.9 or higher ")
endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
endif()

############################################################
# Sanitizer options
############################################################

# Set sanitizers compiler options
if(NOT QE_ENABLE_SANITIZER STREQUAL "none")

  if(QE_ENABLE_SANITIZER STREQUAL "asan" )
    set(GNU_SANITIZER_OPTIONS "-fsanitize=address"
        CACHE STRING "AddressSanitizer Fortran compiler builds." FORCE)
  elseif(QE_ENABLE_SANITIZER STREQUAL "ubsan" )
    set(GNU_SANITIZER_OPTIONS "-fsanitize=undefined"
        CACHE STRING "UndefinedBehaviorSanitizer Fortran compiler builds." FORCE)
  elseif(QE_ENABLE_SANITIZER STREQUAL "tsan" )
    set(GNU_SANITIZER_OPTIONS "-fsanitize=thread"
        CACHE STRING "ThreadSanitizer Fortran compiler builds." FORCE)
  endif()

  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${GNU_SANITIZER_OPTIONS}")
endif()
