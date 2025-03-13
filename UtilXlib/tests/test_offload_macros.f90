program test_gpu_macros
  implicit none
  real a, b
  integer k
#if defined(_OPENACC)
  write(*,*) "Using OpenACC"
  !$acc parallel loop private(a)
#elif defined(__OPENMP_GPU)
  write(*,*) "Using OpenMP GPU offload"
  !$omp target parallel do private(a) map(b)
#elif defined(_OPENMP)
  write(*,*) "Using OpenMP CPU"
#else
  !$omp parallel do private(a)
  write(*,*) "Neither OpenMP nor OpenACC"
#endif
  do k = 1, 1024
    a = b
    a = sin(a)
  enddo
end program
