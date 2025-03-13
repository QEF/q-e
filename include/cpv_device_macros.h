#if defined (_OPENACC)
#define DEV_OMP_NOACC !!!
#define START_WSHARE  !$acc kernels
#define END_WSHARE    !$acc end kernels
#else
#define DEV_OMP_NOACC !$omp
#define START_WSHARE  !$omp workshare
#define END_WSHARE    !$omp end workshare
#endif
