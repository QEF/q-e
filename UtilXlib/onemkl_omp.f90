#if defined(__ONEMKL)

include "mkl_omp_offload.f90"

module onemkl_dfti_gpu
   use mkl_dfti_omp_offload
endmodule onemkl_dfti_gpu

module onemkl_blas_gpu
#if defined(MKL_ILP64)
   use onemkl_blas_omp_offload_ilp64
#else
   use onemkl_blas_omp_offload_lp64
#endif
endmodule onemkl_blas_gpu

module onemkl_blas_no_array_check_gpu
#if defined(MKL_ILP64)
    use onemkl_blas_omp_offload_ilp64_no_array_check
#else
    use onemkl_blas_omp_offload_lp64_no_array_check
#endif
endmodule onemkl_blas_no_array_check_gpu

module onemkl_lapack_gpu
#if defined(MKL_ILP64)
   use onemkl_lapack_omp_offload_ilp64
#else
   use onemkl_lapack_omp_offload_lp64
#endif
endmodule onemkl_lapack_gpu

module onemkl_vsl_gpu
#if defined(MKL_ILP64)
   use onemkl_vsl_omp_offload_ilp64
#else
   use onemkl_vsl_omp_offload_lp64
#endif
endmodule onemkl_vsl_gpu

#elif defined(__OPENMP_GPU)

module onemkl_dfti_gpu
endmodule onemkl_dfti_gpu

module onemkl_blas_gpu
endmodule onemkl_blas_gpu

module onemkl_blas_no_array_check_gpu
endmodule onemkl_blas_no_array_check_gpu

module onemkl_lapack_gpu
endmodule onemkl_lapack_gpu

#endif
