#if defined(__OPENMP_GPU) && defined(__ONEMKL)
include "mkl_dfti_omp_offload.f90"

module onemkl_dfti_omp_offload
   use mkl_dfti_omp_offload
endmodule onemkl_dfti_omp_offload
#endif
