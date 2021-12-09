###########################################################
# Eigensolver_GPU
###########################################################
qe_git_submodule_update(external/eigensolver_gpu)

set(src_eigensolver_gpu
    eigensolver_gpu/lib_eigsolve/dsyevd_gpu.F90
    eigensolver_gpu/lib_eigsolve/dsygst_gpu.F90
    eigensolver_gpu/lib_eigsolve/dsygvdx_gpu.F90
    eigensolver_gpu/lib_eigsolve/dsymv_gpu.F90
    eigensolver_gpu/lib_eigsolve/dsytd2_gpu.F90
    eigensolver_gpu/lib_eigsolve/dsytrd_gpu.F90
    eigensolver_gpu/lib_eigsolve/eigsolve_vars.F90
    eigensolver_gpu/lib_eigsolve/toolbox.F90
    eigensolver_gpu/lib_eigsolve/zheevd_gpu.F90
    eigensolver_gpu/lib_eigsolve/zhegst_gpu.F90
    eigensolver_gpu/lib_eigsolve/zhegvdx_gpu.F90
    eigensolver_gpu/lib_eigsolve/zhemv_gpu.F90
    eigensolver_gpu/lib_eigsolve/zhetd2_gpu.F90
    eigensolver_gpu/lib_eigsolve/zhetrd_gpu.F90)

# See: https://github.com/NVIDIA/Eigensolver_gpu/blob/master/lib_eigsolve/Makefile
# Comment: the flags "-O3 -mp -Mlarge_arrays" are inherited
#          from the global flags of PGI compiler
set(FLAGS -pgf90libs -Mcuda)
set(FLAGS2 -pgf90libs -Mcuda,maxregcount:64)
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(FLAGS3 -pgf90libs -Mcuda=cc35,cc60,nordc,maxregcount:255)
elseif(CMAKE_CUDA_ARCHITECTURES GREATER 60)
    set(FLAGS3 -pgf90libs -Mcuda=cc60,nordc,maxregcount:255)
else()
    set(FLAGS3 -pgf90libs -Mcuda=nordc,maxregcount:255)
endif()
set_source_files_properties(
    eigensolver_gpu/lib_eigsolve/dsyevd_gpu.F90
    eigensolver_gpu/lib_eigsolve/dsygst_gpu.F90
    eigensolver_gpu/lib_eigsolve/dsygvdx_gpu.F90
    eigensolver_gpu/lib_eigsolve/dsytd2_gpu.F90
    eigensolver_gpu/lib_eigsolve/dsytrd_gpu.F90
    eigensolver_gpu/lib_eigsolve/eigsolve_vars.F90
    eigensolver_gpu/lib_eigsolve/toolbox.F90
    eigensolver_gpu/lib_eigsolve/zheevd_gpu.F90
    eigensolver_gpu/lib_eigsolve/zhegst_gpu.F90
    eigensolver_gpu/lib_eigsolve/zhegvdx_gpu.F90
    eigensolver_gpu/lib_eigsolve/zhetrd_gpu.F90
    PROPERTIES COMPILE_OPTIONS "${FLAGS}")
set_source_files_properties(eigensolver_gpu/lib_eigsolve/zhetd2_gpu.F90 eigensolver_gpu/lib_eigsolve/dsymv_gpu.F90
                            PROPERTIES COMPILE_OPTIONS "${FLAGS2}")
set_source_files_properties(eigensolver_gpu/lib_eigsolve/zhemv_gpu.F90 PROPERTIES COMPILE_OPTIONS "${FLAGS3}")

add_library(qe_eigensolver_gpu ${src_eigensolver_gpu})
qe_fix_fortran_modules(qe_eigensolver_gpu)

target_link_libraries(qe_eigensolver_gpu PRIVATE qe_openmp_fortran qe_lapack CUDA::cusolver)

qe_install_targets(qe_eigensolver_gpu)
