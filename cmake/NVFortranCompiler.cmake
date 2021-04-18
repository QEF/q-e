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
    set(QE_CUDA_LINK_OPTIONS ${CUDA_FLAG})

    # CMake default CMAKE_Fortran_FLAGS_RELEASE as -fast -O3
    # -O3 makes the CUDA runs fail at stres_us_gpu.f90, thus override
    set(CMAKE_Fortran_FLAGS_RELEASE "-fast")
endif()
