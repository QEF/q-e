#if defined (_OPENACC) && defined (__CUDA) 
#ifndef __CUF_AND_ACC
#define __CUF_AND_ACC
#endif
#endif 
#if defined (__CUF_AND_ACC) 
#define DEV_ACC !$acc 
#define DEV_OMP !!! 
#define START_WSHARE DEV_ACC  kernels 
#define END_WSHARE   DEV_ACC end  kernels
#else 
#define DEV_ACC !!!
#define DEV_OMP !$omp 
#define START_WSHARE DEV_OMP workshare
#define END_WSHARE   DEV_OMP end workshare
#endif

