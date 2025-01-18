if(QE_ENABLE_OFFLOAD)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM")
        target_compile_options(qe_openmp_fortran INTERFACE "$<$<COMPILE_LANGUAGE:Fortran>:-fopenmp-targets=spir64>")
        target_link_options(qe_openmp_fortran INTERFACE "$<$<LINK_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS};-fopenmp-targets=spir64>")
    else()
        message(FATAL_ERROR "Classic Intel compilers detected. OpenMP offlaod requires Intel OneAPI compilers.")
    endif()
endif()

