#if defined (_OPENACC) 
#ifndef __OPENACC 
#define __OPENACC 
#endif
#endif 
#if defined (__OPENACC) 
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

