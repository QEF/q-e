include(CheckFortranCompilerFlag)

qe_add_global_compile_definitions(__PGI)

# set optimization specific flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mcache_align -Mlarge_arrays")

if(QE_ENABLE_CUDA)
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 20.7)
        set(CUDA_FLAG "-cuda")
    else()
        set(CUDA_FLAG "-Mcuda")
    endif()

    set(QE_CUDA_COMPILE_OPTIONS ${CUDA_FLAG})

    set(GPU_TARGET_COMPILE_OPTIONS)
    if(DEFINED NVFORTRAN_CUDA_VERSION)
        list(APPEND GPU_TARGET_COMPILE_OPTIONS "-gpu=cuda${NVFORTRAN_CUDA_VERSION}")
    endif()

    if(DEFINED NVFORTRAN_CUDA_CC)
        list(APPEND GPU_TARGET_COMPILE_OPTIONS "-gpu=cc${NVFORTRAN_CUDA_CC}")
    endif()

    if(GPU_TARGET_COMPILE_OPTIONS)
      list(APPEND QE_CUDA_COMPILE_OPTIONS ${GPU_TARGET_COMPILE_OPTIONS})
    endif()

    message("   nvfortran CUDA related compile and link options : ${QE_CUDA_COMPILE_OPTIONS}")
    set(CMAKE_REQUIRED_LINK_OPTIONS ${QE_CUDA_COMPILE_OPTIONS})
    check_fortran_compiler_flag("${QE_CUDA_COMPILE_OPTIONS}" NVFORTRAN_CUDA_VALID)
    unset(CMAKE_REQUIRED_LINK_OPTIONS)
    if(NOT NVFORTRAN_CUDA_VALID)
        unset(NVFORTRAN_CUDA_VALID CACHE)
        message(FATAL_ERROR "nvfortran CUDA related option check failed! "
                            "Please check CMakeError.log for the exact error.")
    endif()

    # CMake default CMAKE_Fortran_FLAGS_RELEASE as -fast -O3
    # -O3 makes the CUDA runs fail at stres_us_gpu.f90, thus override
    set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
endif()
