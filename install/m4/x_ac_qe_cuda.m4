# Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
#####
#
# SYNOPSIS
#
# AX_CHECK_CUDA
#
# DESCRIPTION
#
# Simplified compilation for NVidia GPUs using nvhpc compiler
# Assumes a standard installation of a recente nvhpc sdk
# The following variables are substituted in the makefile:
#    cuda_arch
#    cuda_runtime
#    cuda_cflags
#    cuda_fflags
#    cuda_libs
#
# LICENCE
# Public domain
#
#####

AC_DEFUN([X_AC_QE_CUDA], [

# Variables
cuda_arch=
cuda_runtime=
cuda_cflags=
cuda_fflags=
cuda_libs=
   
AC_ARG_WITH([cuda-cc],
   [AS_HELP_STRING([--with-cuda-cc=VAL],[GPU compute capabilities @<:@default=70@:>@])],
   [],
   [with_cuda_cc=70])
   
AC_ARG_WITH([cuda-runtime],
   [AS_HELP_STRING([--with-cuda-runtime=VAL],[CUDA runtime @<:@default=10.2@:>@])],
   [],
   [with_cuda_runtime=10.2])
   
AC_ARG_WITH([cuda-mpi],
   [AS_HELP_STRING([--with-cuda-mpi=VAL],[CUDA-aware MPI (yes|no) @<:@default=no@:>@])],
   [],
   [with_cuda_mpi=no])


AC_ARG_ENABLE([nvtx],
   [AS_HELP_STRING([--enable-nvtx],[Enable compilation for NVTX @<:@default=no@:>@])],
   [enable_nvtx=$enableval],
   [enable_nvtx=no])

if test "x$with_cuda" != "xno"; then
   # NVHPC v.< 21.7 too old (FIXME: still allowing 21.2 for CI)
   if (test "$f90_major_version" -lt 21 ) || (test "$f90_major_version" -eq 21 && test "$f90_minor_version" -lt 2 ) ; then
      AC_MSG_WARN([Compiler version too old, use at least 21.7])
   fi
   if (test "$f90_major_version" -lt 21 ) ; then
      AC_MSG_ERROR([Compiler version too old, use at least 21.7])
   fi
   # NVHPC v. 21.11-22.1 buggy
   if (test "$f90_major_version" -eq 21 && test "$f90_minor_version" -ge 11) ||
      (test "$f90_major_version" -eq 22 && test "$f90_minor_version" -le 1 ) ; then
      AC_MSG_ERROR([Buggy compiler version, upgrade to 22.3 or downgrade to 21.9])
   fi

   # -----------------------------------------
   # Check compiler
   # -----------------------------------------
   AX_CHECK_COMPILE_FLAG([-cuda -gpu=cuda$with_cuda_runtime], [have_cudafor=yes], [have_cudafor=no], [], [MODULE test; use cudafor; END MODULE])
   if test "x$have_cudafor" != "xyes"
   then
      AC_MSG_ERROR([You do not have the cudafor module. Are you using NVHPC compiler?])
   fi
   # -----------------------------------------
   # Headers and libraries
   # -----------------------------------------
   try_dflags="$try_dflags -D__CUDA"
   if test "$use_parallel" -eq 1 && test "$with_cuda_mpi" == "yes"; then 
      try_dflags="$try_dflags -D__GPU_MPI"
   fi
   cuda_libs="-cudalib=cufft,cublas,cusolver,curand"

   cuda_fflags="-cuda -gpu=cc$with_cuda_cc,cuda$with_cuda_runtime"
   # -----------------------------------------
   # Fortran flags
   # -----------------------------------------
   runtime_major_version=`echo $with_cuda_runtime | cut -d. -f1`
   runtime_minor_version=`echo $with_cuda_runtime | cut -d. -f2`
   if test "$runtime_major_version" -lt 10 ||
     (test "$runtime_major_version" -eq 10 && test "$runtime_minor_version" -lt 1 )
   then
       # CUDA toolkit v < 10.1: cusolver not available
       AC_MSG_ERROR([Unsupported CUDA Toolkit, too old])
   fi
   #
   if test "$enable_nvtx" == "yes"; then
      try_dflags="$try_dflags -D__PROFILE_NVTX"
      if test "$runtime_major_version" -gt 12 || \
        (test "$runtime_major_version" -eq 12 && test "$runtime_minor_version" -ge 9); then
         # CUDA >= 12.9: libnvToolsExt.so removed; nvtx3 is header-only
         cuda_libs="${cuda_libs},nvtx"
      else
         # CUDA < 12.9: use legacy shared library
         cuda_fflags="$cuda_fflags -InvToolsExt.h -lnvToolsExt"
      fi
   fi
   # -----------------------------------------
   # C flags 
   # -----------------------------------------   
   cuda_cflags=" -gpu=cc$with_cuda_cc,cuda$with_cuda_runtime"
   ldflags="$ldflags -cuda -gpu=cc$with_cuda_cc,cuda$with_cuda_runtime"
   cuda_arch="$with_cuda_cc"
   cuda_runtime="$with_cuda_runtime"
   ldflags="$ldflags -acc"
   cuda_fflags="$cuda_fflags -acc"
   cuda_cflags="$cuda_cflags -acc"
fi

# Announcing the new variables
AC_SUBST(cuda_arch)
AC_SUBST(cuda_runtime)
AC_SUBST(cuda_cflags)
AC_SUBST(cuda_fflags)
AC_SUBST(cuda_libs)
])
