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
# The following variables are substituted in the makefile:
# NVCC        : the nvcc compiler command.
# NVCCFLAGS   : nvcc specific flags
# CUDA_CFLAGS : CUDA includes
# CUDA_LDLIBS : CUDA libraries
#
# Defines HAVE_CUDA in config.h
#
# LICENCE
# Public domain
#
#####

AC_DEFUN([X_AC_QE_CUDA], [

# Varibles
NVCC=no
CUDA_CFLAGS=
CUDA_LDLIBS=
cuda_fflags=

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
   [AS_HELP_STRING([--with-cuda-runtime=VAL],[CUDA runtime (Pascal: 8.0, Volta: 9.0) @<:@default=8.0@:>@])],
   [],
   [with_cuda_runtime=8.0])



if test "x$with_cuda" != "xno"
then
   # -----------------------------------------
   # Check compiler is PGI
   # -----------------------------------------
   AC_LANG_PUSH([Fortran])
   AC_FC_SRCEXT([f90])
   AX_CHECK_COMPILE_FLAG([-Mcuda], [have_cudafor=yes], [have_cudafor=no], [], [MODULE test; use cudafor; END MODULE])
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
   if test "x$NVCC" = "xno"
   then
      AC_MSG_ERROR([Cannot find nvcc compiler. To enable CUDA, please add path to
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

   if $NVCC -c $NVCCFLAGS conftest.cu &> /dev/null
   then
      ac_compile_nvcc=yes
   fi
   rm -f conftest.cu conftest.o
   AC_MSG_RESULT([$ac_compile_nvcc])

   if test "x$ac_compile_nvcc" = "xno"
   then
      AC_MSG_ERROR([CUDA compiler has problems.])
   fi


   # -----------------------------------------
   # Check for headers and libraries
   # -----------------------------------------
   ax_save_CXXFLAGS="${CXXFLAGS}"
   ax_save_LIBS="${LIBS}"
   

   FCFLAGS="$CUDA_CFLAGS $CXXFLAGS"
   LIBS="$CUDA_LDLIBS $LIBS"

   # And the header and the lib
   AC_CHECK_LIB([cuda], [cuInit], [], AC_MSG_FAILURE([Couldn't find libcuda]))
   AC_CHECK_LIB([cudart], [cudaMalloc], [], AC_MSG_FAILURE([Couldn't find libcudart]))
   AC_CHECK_LIB([cublas], [cublasInit], [], AC_MSG_FAILURE([Couldn't find libcublas]))
   AC_CHECK_LIB([cufft], [cufftPlanMany], [], AC_MSG_FAILURE([Couldn't find libcufft]))
   
   # Returning to the original flags
   CXXFLAGS=${ax_save_CXXFLAGS}
   LIBS=${ax_save_LIBS}

   AC_DEFINE(HAVE_CUDA,1,[Define if we have CUDA])
   try_dflags="$try_dflags -D__CUDA"
   cuda_fflags="-Mcuda=cc$with_cuda_cc,cuda$with_cuda_runtime"
   ldflags="$ldflags -Mcuda=cc$with_cuda_cc,cuda$with_cuda_runtime"
   gpu_arch="$with_cuda_cc"
   gpu_runtime="$with_cuda_runtime"
   pgi_cuda_libs="-Mcudalib=cufft,cublas,cusolver \$(TOPDIR)/EIGENSOLVER_GPU/lib_eigsolve/lib_eigsolve.a"
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
AC_SUBST(pgi_cuda_libs)
])
