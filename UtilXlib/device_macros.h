
#if defined(_CUDA)
#  ifndef __CUDA
#     define __CUDA
#  endif
#endif

#if defined(__CUDA) || defined(__OPENACC) || defined(__OPENMP5)
#  define __HAVE_DEVICE
#endif

