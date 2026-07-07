#--------------------------------------------------------------------
# QE_GPU option
#--------------------------------------------------------------------

# Guess QE_GPU_DEFAULT if QE_GPU_ARCHS is provided. Only use the first entry if there are more.
if(QE_GPU_ARCHS)
  list(GET QE_GPU_ARCHS 0 GPU_ARCH_0)
  if(GPU_ARCH_0 MATCHES "^sm_")
    set(QE_GPU_DEFAULT "openacc;cuda")
  elseif(GPU_ARCH_0 MATCHES "^gfx")
    set(QE_GPU_DEFAULT "openmp;rocm")
  elseif(GPU_ARCH_0 MATCHES "^intel_gpu_")
    set(QE_GPU_DEFAULT "openmp;oneapi")
  else()
    set(QE_GPU_DEFAULT "")
  endif()
endif()

set(VALID_QE_GPU_FEATURES "openacc" "openmp" "cuda" "rocm" "oneapi")
set(QE_GPU ${QE_GPU_DEFAULT} CACHE STRING "Semicolon-separated list of GPU features to enable (openmp,cuda,rocm,oneapi)")
set_property(CACHE QE_GPU PROPERTY STRINGS ${VALID_QE_GPU_FEATURES})

set(QE_ENABLE_OPENACC OFF)
set(QE_ENABLE_OFFLOAD OFF)
set(QE_ENABLE_CUDA OFF)
set(QE_ENABLE_ROCM OFF)
set(QE_ENABLE_ONEAPI OFF)

# Perform QE_GPU option check
foreach(GPU_FEATURE IN LISTS QE_GPU)
  # verify the entry
  if(NOT GPU_FEATURE IN_LIST VALID_QE_GPU_FEATURES)
    message(FATAL_ERROR "Invalid QE_GPU selection \"${GPU_FEATURE}\", valid values are \"${VALID_QE_GPU_FEATURES}\"")
  endif()

  if(GPU_FEATURE STREQUAL "openacc")
    set(QE_ENABLE_OPENACC ON)
  endif()
  if(GPU_FEATURE STREQUAL "openmp")
    set(QE_ENABLE_OFFLOAD ON)
  endif()
  if(GPU_FEATURE STREQUAL "cuda")
    set(QE_ENABLE_CUDA ON)
  endif()
  if(GPU_FEATURE STREQUAL "rocm")
    set(QE_ENABLE_ROCM ON)
  endif()
  if(GPU_FEATURE STREQUAL "oneapi")
    set(QE_ENABLE_ONEAPI ON)
  endif()
endforeach()

if("cuda" IN_LIST QE_GPU AND "rocm" IN_LIST QE_GPU OR
   "cuda" IN_LIST QE_GPU AND "oneapi" IN_LIST QE_GPU OR
   "rocm" IN_LIST QE_GPU AND "oneapi" IN_LIST QE_GPU)
    message(FATAL_ERROR "Invalid QE_GPU selection\"${QE_GPU}\", \"cuda\", \"rocm\" and \"oneapi\" don't work together. Select one!")
endif()

if(QE_GPU)
  message(STATUS "Enable GPU features QE_GPU=${QE_GPU}")
endif()
