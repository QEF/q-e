qe_add_global_compile_definitions(__PGI)

# set optimization specific flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mcache_align -Mlarge_arrays")

if(QE_ENABLE_CUDA)
    set(QE_CUDA_COMPILE_OPTIONS "-Mcuda")
    set(QE_CUDA_LINK_OPTIONS "-Mcuda")
endif()