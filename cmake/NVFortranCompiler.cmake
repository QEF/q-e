include(CheckFortranCompilerFlag)

qe_add_global_compile_definitions(__PGI)

# set optimization specific flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mcache_align -Mlarge_arrays -Mbackslash")

if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 21.11 AND CMAKE_Fortran_COMPILER_VERSION VERSION_LESS_EQUAL 22.2)
    if(QE_ENABLE_OPENACC AND QE_ENABLE_OPENMP)
        message(FATAL_ERROR "NVHPC 21.11-22.2 have a severe bug causing hanging in runs"
                            " when QE is compiled with both OpenMP and OpenACC. "
                            "Use a different compiler release or dislable OpenMP with potential performance loss.")
    endif()
endif()

# set up GPU architecture options which can be applied to CUDA Fortran, OpenACC and OpenMP offload
set(GPU_TARGET_COMPILE_OPTIONS)
if(QE_ENABLE_CUDA OR QE_ENABLE_OPENACC OR QE_ENABLE_OFFLOAD)
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS 19.10)
        message(FATAL_ERROR "Compiler Version ${CMAKE_Fortran_COMPILER_VERSION}. "
                            "GPU acceleration requires PGI 19.10 or NVIDIA HPC SDK 20.7 or higher!")
    endif()

    if(DEFINED NVFORTRAN_CUDA_VERSION)
        list(APPEND GPU_TARGET_COMPILE_OPTIONS "-gpu=cuda${NVFORTRAN_CUDA_VERSION}")
    endif()

    if(DEFINED NVFORTRAN_CUDA_CC)
        list(APPEND GPU_TARGET_COMPILE_OPTIONS "-gpu=cc${NVFORTRAN_CUDA_CC}")
    elseif(DEFINED QE_GPU_ARCHS)
        string(REPLACE "sm_" "" CUDA_ARCH_NUMBERS "${QE_GPU_ARCHS}")
        string(REPLACE ";" ",cc" OFFLOAD_ARCH "${CUDA_ARCH_NUMBERS}")
        list(APPEND GPU_TARGET_COMPILE_OPTIONS "-gpu=cc${OFFLOAD_ARCH}")
    endif()
endif()

if(QE_ENABLE_CUDA)
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 20.7)
        set(CUDA_FLAG "-cuda")
    else()
        set(CUDA_FLAG "-Mcuda")
    endif()

    set(QE_CUDA_COMPILE_OPTIONS ${CUDA_FLAG})

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

if(QE_ENABLE_OPENACC)
    if(GPU_TARGET_COMPILE_OPTIONS)
        target_compile_options(qe_openacc_fortran INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:${GPU_TARGET_COMPILE_OPTIONS}>")
        target_compile_options(qe_openacc_c INTERFACE "$<$<COMPILE_LANGUAGE:C>:${GPU_TARGET_COMPILE_OPTIONS}>")
    endif()
endif()

if(QE_ENABLE_OFFLOAD)
    target_compile_options(qe_openmp_fortran INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-mp=gpu>")
    target_link_options(qe_openmp_fortran INTERFACE "$<$<LINK_LANGUAGE:Fortran>:-mp=gpu>")
    if(GPU_TARGET_COMPILE_OPTIONS)
        target_compile_options(qe_openmp_fortran INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:${GPU_TARGET_COMPILE_OPTIONS}>")
    endif()
endif()
