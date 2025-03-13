###########################################################
# qe_fftx target
###########################################################
set(f_src_fftx
    fft_scatter.f90
    fft_scatter_2d.f90
    scatter_mod.f90 
    fft_ggen.f90
    fft_fwinv.f90
    fft_scalar.f90
    fftw_interfaces.f90
    fft_parallel.f90
    fft_parallel_2d.f90
    fft_interfaces.f90
    fft_interpolate.f90
    stick_base.f90
    fft_smallbox.f90
    fft_smallbox_type.f90
    fft_support.f90
    fft_error.f90
    fft_types.f90
    tg_gather.f90
    fft_helper_subroutines.f90
    fft_param.f90)
if(QE_ENABLE_CUDA)
    set(f_src_fftx
        ${f_src_fftx}
        fft_scalar.cuFFT.f90
        fft_buffers.f90
        fft_scatter_gpu.f90
        fft_scatter_2d_gpu.f90)
endif()
qe_enable_cuda_fortran("${f_src_fftx}")

set(c_src_fftx
    fft_stick.c
    fftw.c
    fftw_sp.c
    fftw_dp.c)

add_library(qe_fftx_c ${c_src_fftx})

qe_add_library(qe_fftx 
    ${f_src_fftx} 
    ${qe_ext_fft_wrappers})

target_link_libraries(qe_fftx
    PRIVATE
        qe_fftx_c
        qe_ext_fft
        qe_openmp_fortran
        qe_mpi_fortran)
if(QE_ENABLE_CUDA)
    target_link_libraries(qe_fftx
        PUBLIC
            qe_openacc_fortran
        PRIVATE
            CUDA::cufft)
endif()
        
qe_install_targets(qe_fftx qe_fftx_c)
