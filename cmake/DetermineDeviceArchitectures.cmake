# derive AMD GPU architectures from CMAKE_HIP_ARCHITECTURES or detect them using amdgpu-arch (LLVM15+)
function(detectAMDGPU)
  if(CMAKE_HIP_ARCHITECTURES)
    set(QE_GPU_ARCHS_DETECTED ${CMAKE_HIP_ARCHITECTURES})
  else()
    find_program(AMDGPU_ARCH_EXE amdgpu-arch)
    if(AMDGPU_ARCH_EXE)
      execute_process(COMMAND ${AMDGPU_ARCH_EXE} OUTPUT_VARIABLE AMD_GPU_ARCH)
      string(REGEX REPLACE "\n" ";" AMD_GPU_ARCH "${AMD_GPU_ARCH}")
      if(AMD_GPU_ARCH)
        list(APPEND QE_GPU_ARCHS_DETECTED ${AMD_GPU_ARCH})
      else()
        message("amdgpu-arch didn't find AMD GPUs. Cannot auto-detect AMD GPU architectures.")
      endif()
    else()
      message("amdgpu-arch not avaible. Cannot auto-detect AMD GPU architectures.")
    endif()
  endif()
  set(QE_GPU_ARCHS
      ${QE_GPU_ARCHS_DETECTED}
      PARENT_SCOPE)
endfunction()

# check QE_GPU_ARCHS and CMAKE_HIP_ARCHITECTURES consistency
function(verifyAMDGPUconsistency)
  if(CMAKE_HIP_ARCHITECTURES AND NOT CMAKE_HIP_ARCHITECTURES STREQUAL QE_GPU_ARCHS)
    message(FATAL_ERROR "CMAKE_HIP_ARCHITECTURES=${CMAKE_HIP_ARCHITECTURES} doesn't match ${QE_GPU_ARCHS}.")
  endif()
endfunction()

# derive NVIDIA GPU architectures from CMAKE_CUDA_ARCHITECTURES or detect them using nvptx-arch (LLVM16+)
function(detectNVIDIAGPU)
  if(CMAKE_CUDA_ARCHITECTURES)
    set(QE_GPU_ARCHS_DETECTED)
    foreach(CUDA_ARCH_NUM IN LISTS CMAKE_CUDA_ARCHITECTURES)
      list(APPEND QE_GPU_ARCHS_DETECTED "sm_${CUDA_ARCH_NUM}")
    endforeach()
  else()
    find_program(NVPTX_ARCH_EXE nvptx-arch)
    if(NVPTX_ARCH_EXE)
      execute_process(COMMAND ${NVPTX_ARCH_EXE} OUTPUT_VARIABLE NVIDIA_GPU_ARCH)
      string(REGEX REPLACE "\n" ";" NVIDIA_GPU_ARCH "${NVIDIA_GPU_ARCH}")
      if(NVIDIA_GPU_ARCH)
        list(APPEND QE_GPU_ARCHS_DETECTED ${NVIDIA_GPU_ARCH})
      else()
        message("nvptx-arch didn't find NVIDIA GPUs. Cannot auto-detect NVIDIA GPU architectures.")
      endif()
    else()
      message("nvptx-arch not avaible. Cannot auto-detect NVIDIA GPU architectures.")
    endif()
  endif()
  set(QE_GPU_ARCHS
      ${QE_GPU_ARCHS_DETECTED}
      PARENT_SCOPE)
endfunction()

# check QE_GPU_ARCHS and CMAKE_CUDA_ARCHITECTURES consistency
function(verifyNVIDIAGPUconsistency)
  string(REPLACE "sm_" "" CUDA_ARCH_NUMBERS "${QE_GPU_ARCHS}")
  if(CMAKE_CUDA_ARCHITECTURES AND NOT CMAKE_CUDA_ARCHITECTURES STREQUAL CUDA_ARCH_NUMBERS)
    message(FATAL_ERROR "CMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES} doesn't match ${QE_GPU_ARCHS}.")
  endif()
endfunction()

# auto detect QE_GPU_ARCHS if not set by user and GPU features are enabled.
# CMAKE_CUDA/HIP_ARCHITECTURES are used as hints
if(NOT QE_GPU_ARCHS AND (QE_ENABLE_CUDA OR ENABLE_ROCM))
  if(QE_ENABLE_CUDA)
    detectNVIDIAGPU()
  endif()
  if(QE_ENABLE_ROCM)
    detectAMDGPU()
  endif()

  if(NOT QE_GPU_ARCHS)
    message(
      WARNING
        "QE_GPU_ARCHS was neither set nor derivable and auto-detection failed. "
	"Some compilers may decide to target all the architectures they support and "
	"hence slow down the overall compilation process. ")
  endif()
endif()

list(REMOVE_DUPLICATES QE_GPU_ARCHS)

# make sure QE_GPU_ARCHS is consistent with CMAKE_HIP_ARCHITECTURES or CMAKE_CUDA_ARCHITECTURES.
if(QE_ENABLE_CUDA)
  verifyNVIDIAGPUconsistency()
endif()
if(QE_ENABLE_ROCM)
    verifyAMDGPUconsistency()
endif()

set(QE_GPU_ARCHS
    ${QE_GPU_ARCHS}
    CACHE STRING "Accelerator device architectures" FORCE)

if(QE_GPU_ARCHS)
  message(STATUS "GPU device architectures: ${QE_GPU_ARCHS}")
endif()

# QE_GPU_ARCHS is the single source of truth and thus overwrite CMAKE_CUDA/HIP_ARCHITECTURES
if(QE_ENABLE_CUDA)
  string(REPLACE "sm_" "" CUDA_ARCH_NUMBERS "${QE_GPU_ARCHS}")
  set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_NUMBERS} CACHE STRING "CUDA architectures" FORCE)
endif()
if(QE_ENABLE_ROCM)
  set(CMAKE_HIP_ARCHITECTURES ${QE_GPU_ARCHS} CACHE STRING "HIP architectures" FORCE)
endif()
