program test_gpu_macros
  implicit none
  real a
  integer k
#if defined(_OPENACC)
  write(*,*) "Using OpenACC"
  !$acc parallel loop private(a)
#elif defined(__OPENMP_GPU)
  write(*,*) "Using OpenMP GPU offload"
  !$omp target teams distribute parallel do private(a)
#elif defined(_OPENMP)
  write(*,*) "Using OpenMP CPU"
#else
  !$omp parallel do private(a)
  write(*,*) "Neither OpenMP nor OpenACC"
#endif
  do k = 1, 1024
    a = k
    a = sin(a)
  enddo
end program
