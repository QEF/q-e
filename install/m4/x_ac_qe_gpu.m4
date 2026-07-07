# Copyright (C) 2026 Quantum ESPRESSO Foundation
#

AC_DEFUN([X_AC_QE_GPU], [

AC_ARG_WITH(gpu,
   [AS_HELP_STRING([--with-gpu],
       [(cuda|openmp) Use "cuda" for NVIDIA, "omp5" for AMD or INTEL GPUs (default:cuda)])],
   [if test "$withval" = "openmp" ; then
      with_gpu=2
   else
      with_gpu=1
   fi],
   [with_gpu=0])
 
AC_ARG_WITH(cuda,
   [AS_HELP_STRING([--with-cuda],
       [obsolete, use --with-gpu=cuda instead])],
   [if test "$withval" != "no" ; then
       AC_MSG_WARN([--with-cuda is obsolete, use --with-gpu=cuda instead])
       with_gpu=1
   fi],
   [])
# do not set the default with_gpu=0 here, it breaks the result of --with-gpu 
   
])
