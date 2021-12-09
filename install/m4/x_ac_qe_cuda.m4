# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation
#####
#
# SYNOPSIS
#
# AX_CHECK_CUDA
#
# DESCRIPTION
#
# Figures out if CUDA Driver API/nvcc is available, i.e. existence of:
#   nvcc
#   cuda.h
#   libcuda.a
#
# If something isn't found, fails straight away.
#
# The following variables (for C) are substituted in the makefile:
#    NVCC        : the nvcc compiler command.
#    NVCCFLAGS   : nvcc specific flags
#    CUDA_CFLAGS : CUDA includes
#    CUDA_LDLIBS : CUDA libraries
# For Fortran:
#    gpu_arch
#    gpu_runtime
#    cuda_fflags
#    cuda_libs
#    cuda_extlibs
#    cuda_path
#
# LICENCE
# Public domain
#
#####

AC_DEFUN([X_AC_QE_CUDA], [

# Variables
NVCC=no
CUDA_CFLAGS=
CUDA_LDLIBS=
cuda_path=
cuda_fflags=
cuda_libs=
cuda_extlibs=

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

AC_ARG_ENABLE([cuda-env-check],
   [AS_HELP_STRING([--enable-cuda-env-check],[The configure script will check CUDA installation and report problems @<:@default=no@:>@])],
   [],
   [enable_cuda_env_check=no])

if test "$f90_major_version" -gt 20 || (test "$f90_major_version" -eq 20 && test "$f90_minor_version" -ge 7); then
   mMcuda="-cuda -gpu"
   mMcudalib="-cudalib"
else
   mMcuda="-Mcuda"
   mMcudalib="-Mcudalib"
fi

if test "x$with_cuda" != "xno"
then
   # -----------------------------------------
   # Check compiler is PGI
   # -----------------------------------------
   AC_LANG_PUSH([Fortran])
   AC_FC_SRCEXT([f90])
   AX_CHECK_COMPILE_FLAG([$mMcuda=cuda$with_cuda_runtime], [have_cudafor=yes], [have_cudafor=no], [], [MODULE test; use cudafor; END MODULE])
   AC_LANG_POP([Fortran])
   if test "x$have_cudafor" != "xyes"
   then
      AC_MSG_ERROR([You do not have the cudafor module. Are you using a PGI compiler?])
   fi
   # -----------------------------------------
   # Setup CUDA paths
   # -----------------------------------------
   if test "x$with_cuda" != "xyes"
   then
      #AX_NORMALIZE_PATH([with_cuda], ["/"])
      CUDAPATH="$with_cuda"
      CUDA_CFLAGS+=" -I$with_cuda/include"
      CUDA_LDLIBS+=" -L$with_cuda/lib64"
   else
      AC_CHECK_FILE(/usr/local/cuda/,[CUDAPATH="/usr/local/cuda"],[])
      AC_CHECK_FILE(/usr/local/cuda/include,[CUDA_CFLAGS+=" -I/usr/local/cuda/include"],[CUDA_CFLAGS=""])
      AC_CHECK_FILE(/usr/local/cuda/lib64,[CUDA_LDLIBS+=" -L/usr/local/cuda/lib64"],[])
   fi
   CUDA_LDLIBS+=" -lcuda -lcudart -lcublas -lcufft"



   # -----------------------------------------
   # Checking for nvcc
   # -----------------------------------------
   AC_PATH_PROG([NVCC],[nvcc],[no],[$PATH:$CUDAPATH/bin])
   if test "x$NVCC" = "xno" && test "x$enable_cuda_env_check" = "xyes"
   then
      AC_MSG_WARN([Cannot find nvcc compiler. To enable CUDA, please add path to
                    nvcc in the PATH environment variable and/or specify the path
                    where CUDA is installed using: --with-cuda=PATH])
   fi


   # -----------------------------------------
   # Setup nvcc flags
   # -----------------------------------------
   AC_ARG_VAR(NVCCFLAGS,[Additional nvcc flags (example: NVCCFLAGS="-arch=compute_30 -code=sm_30")])
   if test x$DEBUG = xtrue
   then
      NVCCFLAGS+=" -g -arch=compute_$with_cuda_cc -code=sm_$with_cuda_cc"
   else
      NVCCFLAGS+=" -O3 -arch=compute_$with_cuda_cc -code=sm_$with_cuda_cc"
   fi
   

   # -----------------------------------------
   # Check if nvcc works
   # -----------------------------------------
   ac_compile_nvcc=no
   AC_MSG_CHECKING([whether nvcc works])
   cat>conftest.cu<<EOF
   __global__ static void test_cuda() {
      const int tid = threadIdx.x;
      const int bid = blockIdx.x;
      __syncthreads();
   }
EOF

   if test "x$NVCC" != "xno"
   then
      if $NVCC -c $NVCCFLAGS conftest.cu &> /dev/null
      then
         ac_compile_nvcc=yes
      fi
   fi
   rm -f conftest.cu conftest.o
   AC_MSG_RESULT([$ac_compile_nvcc])

   if test "x$ac_compile_nvcc" = "xno" && test "x$enable_cuda_env_check" = "xyes"
   then
      AC_MSG_WARN([CUDA compiler has problems.])
   fi


   # -----------------------------------------
   # Check for headers and libraries
   # -----------------------------------------
   ax_save_CXXFLAGS="${CXXFLAGS}"
   ax_save_LIBS="${LIBS}"
   

   FCFLAGS="$CUDA_CFLAGS $CXXFLAGS"
   LIBS="$CUDA_LDLIBS $LIBS"

   # And the header and the lib
   if test "x$enable_cuda_env_check" = "xyes"
   then
      AC_CHECK_LIB([cuda], [cuInit], [], AC_MSG_FAILURE([Couldn't find libcuda]))
      AC_CHECK_LIB([cudart], [cudaMalloc], [], AC_MSG_FAILURE([Couldn't find libcudart]))
      AC_CHECK_LIB([cublas], [cublasInit], [], AC_MSG_FAILURE([Couldn't find libcublas]))
      AC_CHECK_LIB([cufft], [cufftPlanMany], [], AC_MSG_FAILURE([Couldn't find libcufft]))
   fi

   new_cusolver="yes"
   CC_stash=$CC
   CC=nvcc
   AC_CHECK_LIB([cusolver], [cusolverDnZhegvdx_bufferSize], [], new_cusolver="no")
   CC=$CC_stash
   
   # Returning to the original flags
   CXXFLAGS=${ax_save_CXXFLAGS}
   LIBS=${ax_save_LIBS}

   try_dflags="$try_dflags -D__CUDA"
   cuda_extlibs="devxlib"
   cuda_libs="$mMcudalib=cufft,cublas,cusolver,curand \$(TOPDIR)/external/devxlib/src/libdevXlib.a"
   cuda_fflags="$mMcuda=cc$with_cuda_cc,cuda$with_cuda_runtime"
   cuda_fflags="$cuda_fflags \$(MOD_FLAG)\$(TOPDIR)/external/devxlib/src"
   cuda_fflags="$cuda_fflags \$(MOD_FLAG)\$(TOPDIR)/external/devxlib/include"
   
   if test ${new_cusolver} != yes; then
       cuda_fflags="$cuda_fflags \$(MOD_FLAG)\$(TOPDIR)/EIGENSOLVER_GPU/lib_eigsolve"
       cuda_extlibs="$cuda_extlibs eigensolver"
       cuda_libs="$cuda_libs \$(TOPDIR)/EIGENSOLVER_GPU/lib_eigsolve/lib_eigsolve.a"
       AC_MSG_WARN([Using legacy custom solver.])
   else
       try_dflags="$try_dflags -D__USE_CUSOLVER"
   fi
   ldflags="$ldflags $mMcuda=cc$with_cuda_cc,cuda$with_cuda_runtime"
   gpu_arch="$with_cuda_cc"
   gpu_runtime="$with_cuda_runtime"
   cuda_path="$CUDAPATH"
   if test "$enable_openacc" == "yes"; then
      ldflags="$ldflags -acc"
      cuda_fflags="$cuda_fflags -acc"
      CUDA_CFLAGS="$CUDA_CFLAGS -acc -gpu=cc$with_cuda_cc,cuda$with_cuda_runtime"
   fi

fi

# Announcing the new variables
# For C (maybe needed in the future)
AC_SUBST([NVCC])
AC_SUBST([NVCCFLAGS])
AC_SUBST([CUDA_CFLAGS])
AC_SUBST([CUDA_LDLIBS])
# And for Fortran
AC_SUBST(gpu_arch)
AC_SUBST(gpu_runtime)
AC_SUBST(cuda_fflags)
AC_SUBST(cuda_libs)
AC_SUBST(cuda_extlibs)
AC_SUBST(cuda_path)
])
