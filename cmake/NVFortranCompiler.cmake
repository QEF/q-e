include(CheckFortranSourceCompiles)

qe_add_global_compile_definitions(__PGI)

# set optimization specific flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mcache_align -Mlarge_arrays")

if(QE_ENABLE_OPENACC)
    add_library(OpenACC::OpenACC_Fortran INTERFACE IMPORTED)
    set_target_properties(OpenACC::OpenACC_Fortran PROPERTIES
                          INTERFACE_COMPILE_OPTIONS "-acc"
                          INTERFACE_LINK_OPTIONS "-acc")
endif()

if(QE_ENABLE_CUDA)
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 20.7)
        set(CUDA_FLAG "-cuda")
    else()
        set(CUDA_FLAG "-Mcuda")
    endif()

    set(QE_CUDA_COMPILE_OPTIONS ${CUDA_FLAG})
    set(QE_CUDA_LINK_OPTIONS ${CUDA_FLAG})

    if(DEFINED NVFORTRAN_CUDA_VERSION)
        set(CMAKE_REQUIRED_FLAGS "-Mfree ${CUDA_FLAG} -gpu=cuda${NVFORTRAN_CUDA_VERSION}")
        check_fortran_source_compiles("program abc; end program" NVFORTRAN_CUDA_VERSION_VALID)
        unset(CMAKE_REQUIRED_FLAGS)
        if(NOT NVFORTRAN_CUDA_VERSION_VALID)
            unset(NVFORTRAN_CUDA_VERSION_VALID CACHE)
            message(FATAL_ERROR "nvfortran CUDA version check failed! "
                                "NVFORTRAN_CUDA_VERSION=${NVFORTRAN_CUDA_VERSION} (-gpu=cuda${NVFORTRAN_CUDA_VERSION}) not accepted")
        endif()
        list(APPEND QE_CUDA_COMPILE_OPTIONS "-gpu=cuda${NVFORTRAN_CUDA_VERSION}")
    endif()

    if(DEFINED NVFORTRAN_CUDA_CC)
        set(CMAKE_REQUIRED_FLAGS "-Mfree ${CUDA_FLAG} -gpu=cc${NVFORTRAN_CUDA_CC}")
        check_fortran_source_compiles("program abc; end program" NVFORTRAN_CUDA_CC_VALID)
        unset(CMAKE_REQUIRED_FLAGS)
        if(NOT NVFORTRAN_CUDA_CC_VALID)
            unset(NVFORTRAN_CUDA_CC_VALID CACHE)
            message(FATAL_ERROR "nvfortran GPU architecture check failed! "
                                "NVFORTRAN_CUDA_CC=${NVFORTRAN_CUDA_CC} (-gpu=cc${NVFORTRAN_CUDA_CC}) not accepted")
        endif()
        list(APPEND QE_CUDA_COMPILE_OPTIONS "-gpu=cc${NVFORTRAN_CUDA_CC}")
    endif()

    message("   nvfortran CUDA related compile and link options : ${QE_CUDA_COMPILE_OPTIONS}")
    # CMake default CMAKE_Fortran_FLAGS_RELEASE as -fast -O3
    # -O3 makes the CUDA runs fail at stres_us_gpu.f90, thus override
    set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
endif()
