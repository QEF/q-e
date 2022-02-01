# Copyright (C) 2001-2022 Quantum ESPRESSO Foundation
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
#    gpu_arch
#    cuda_runtime
#    cuda_cflags
#    cuda_fflags
#    cuda_libs
#    cuda_extlibs
#
# LICENCE
# Public domain
#
#####

AC_DEFUN([X_AC_QE_CUDA], [

# Variables
gpu_arch=
cuda_runtime=
cuda_cflags=
cuda_fflags=
cuda_libs=
# FIXME: currently devxlib is needed also for non-CUDA compilation
cuda_extlibs=devxlib

# Provide your CUDA path with this
AC_ARG_WITH([cuda],
   [AS_HELP_STRING([--with-cuda=PATH],[prefix where CUDA is installed @<:@default=no@:>@])],
   [],
   [with_cuda=no])
   
AC_ARG_WITH([cuda-cc],
   [AS_HELP_STRING([--with-cuda-cc=VAL],[GPU architecture (Kepler: 35, Pascal: 60, Volta: 70) @<:@default=35@:>@])],
   [],
   [with_cuda_cc=35])
   
AC_ARG_WITH([cuda-runtime],
   [AS_HELP_STRING([--with-cuda-runtime=VAL],[CUDA runtime (Pascal: 8+, Volta: 9+) @<:@default=10.1@:>@])],
   [],
   [with_cuda_runtime=10.1])

AC_ARG_ENABLE([openacc],
   [AS_HELP_STRING([--enable-openacc],[Enable compilation with OPENACC @<:@default=yes@:>@])],
   [],
   [enable_openacc=yes])

if test "$f90_major_version" -gt 20 || (test "$f90_major_version" -eq 20 && test "$f90_minor_version" -ge 7); then
   # NVHPC v. 20.7 and later
   mMcuda="-cuda -gpu"
   mMcudalib="-cudalib"
else
   # NVHPC previous to v. 20.7
   mMcuda="-Mcuda"
   mMcudalib="-Mcudalib"
fi

if test "x$with_cuda" != "xno"
then
   # -----------------------------------------
   # Check compiler
   # -----------------------------------------
   AC_LANG_PUSH([Fortran])
   AC_FC_SRCEXT([f90])
   AX_CHECK_COMPILE_FLAG([$mMcuda=cuda$with_cuda_runtime], [have_cudafor=yes], [have_cudafor=no], [], [MODULE test; use cudafor; END MODULE])
   AC_LANG_POP([Fortran])
   if test "x$have_cudafor" != "xyes"
   then
      AC_MSG_ERROR([You do not have the cudafor module. Are you using NVHPC compiler?])
   fi
   # -----------------------------------------
   # Headers and libraries
   # -----------------------------------------
   try_dflags="$try_dflags -D__CUDA"
   cuda_extlibs="devxlib"
   cuda_libs="$mMcudalib=cufft,cublas,cusolver,curand \$(TOPDIR)/external/devxlib/src/libdevXlib.a"
   
   cuda_fflags="$mMcuda=cc$with_cuda_cc,cuda$with_cuda_runtime"
   cuda_fflags="$cuda_fflags \$(MOD_FLAG)\$(TOPDIR)/external/devxlib/src"
   cuda_fflags="$cuda_fflags \$(MOD_FLAG)\$(TOPDIR)/external/devxlib/include"
   # -----------------------------------------
   # Fortran flags
   # -----------------------------------------   
   runtime_major_version=`echo $with_cuda_runtime | cut -d. -f1`
   runtime_minor_version=`echo $with_cuda_runtime | cut -d. -f2`
   if test "$runtime_major_version" -lt 10 || 
         ( "$runtime_major_version" -eq 10 && "$runtime_minor_version" -lt 1 )
   then
       # CUDA toolkit v < 10.1: new solver not available
       cuda_fflags="$cuda_fflags \$(MOD_FLAG)\$(TOPDIR)/EIGENSOLVER_GPU/lib_eigsolve"
       cuda_extlibs="$cuda_extlibs eigensolver"
       cuda_libs="$cuda_libs \$(TOPDIR)/EIGENSOLVER_GPU/lib_eigsolve/lib_eigsolve.a"
       AC_MSG_WARN([Using legacy custom solver.])
   else
       try_dflags="$try_dflags -D__USE_CUSOLVER"
   fi
   # -----------------------------------------
   # C flags - not sure whether they are suitable for old version as well
   # -----------------------------------------   
   cuda_cflags=" -I$with_cuda/include -gpu=cc$with_cuda_cc,cuda$with_cuda_runtime"
   ldflags="$ldflags $mMcuda=cc$with_cuda_cc,cuda$with_cuda_runtime"
   gpu_arch="$with_cuda_cc"
   cuda_runtime="$with_cuda_runtime"
   if test "$enable_openacc" == "yes"; then
      ldflags="$ldflags -acc"
      cuda_fflags="$cuda_fflags -acc"
      cuda_cflags="$cuda_cflags -acc"
   fi

fi

# Announcing the new variables
AC_SUBST(gpu_arch)
AC_SUBST(cuda_runtime)
AC_SUBST(cuda_cflags)
AC_SUBST(cuda_fflags)
AC_SUBST(cuda_libs)
AC_SUBST(cuda_extlibs)
])
