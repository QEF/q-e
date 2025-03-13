# Check compiler version
if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 4.9)
  message(FATAL_ERROR "Requires GCC 4.9 or higher ")
endif()

if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
endif()

if(QE_ENABLE_OFFLOAD)
  if(NOT DEFINED QE_GPU_ARCHS)
    message(FATAL_ERROR "Requires QE_GPU_ARCHS option. For example, sm_80 for NVIDIA A100 or gfx90a for AMD MI250X.")
  endif()

  if(QE_GPU_ARCHS MATCHES "sm_")
    set(OFFLOAD_TARGET nvptx-none)
  elseif(QE_GPU_ARCHS MATCHES "gfx")
    set(OFFLOAD_TARGET amdgcn-amdhsa)
  else()
    message(FATAL_ERROR "Cannot derive OFFLOAD_TARGET from QE_GPU_ARCHS.")
  endif()
  target_compile_options(qe_openmp_fortran INTERFACE "-foffload=${OFFLOAD_TARGET};-foffload-options=-lm -latomic")

  if(OFFLOAD_TARGET STREQUAL "amdgcn-amdhsa")
    target_compile_options(qe_openmp_fortran INTERFACE "-foffload-options=${OFFLOAD_TARGET}=-march=${QE_GPU_ARCHS}")
  elseif(OFFLOAD_TARGET STREQUAL "nvptx-none")
    target_compile_options(qe_openmp_fortran INTERFACE "-foffload-options=${OFFLOAD_TARGET}=-misa=${QE_GPU_ARCHS}")
  endif()

  target_link_options(qe_openmp_fortran INTERFACE "$<$<LINK_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS}>")
else()
  target_compile_options(qe_openmp_fortran INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-foffload=disable>")
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
