  !
  ! Copyright (C) 2023-2026 EPW-Collaboration
#if defined (__ONEMKL)
  include "mkl_blas_omp_offload.f90"
#endif

MODULE epw_onemkl
#if defined (__ONEMKL)
  USE onemkl_blas_omp_offload_lp64
#endif
  implicit none
END MODULE epw_onemkl
