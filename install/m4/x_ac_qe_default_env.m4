# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_DEFAULT_ENV], [

# Non-standard precious variables
AC_ARG_VAR(EXTLIB_FLAGS, This variable controls the flags passed to internal BLAS and LAPACK libraries)

# store variables from the environment, if set (may be or not be set)
# If set, they take precedence over configure internal choice.
# Flags and libraries are accepted without further testing;
# compilers are tested. Specify compiler name only, not the full path
# (i.e. F90=/usr/local/bin/f90 may not work, use F90=f90)

topdir=$TOPDIR     # current directory
arch=$ARCH         # see below for recognized architectures
env_cc=$CC         # C compiler (must be in the execution path)
cpp=$CPP           # C preprocessor (as above)
cflags=$CFLAGS     # Flags for C compiler
cppflags=$CPPFLAGS # Flags for C preprocessor
dflags=$DFLAGS     # Fortran file preprocessing options, e.g. -D__DEFINE_THIS
iflags=$IFLAGS     # Location of include files - shouldn't be needed
f77=$F77           # Fortran 77 serial compiler (must be in execution path)
f90=$F90           # Fortran 90 serial compiler (must be in execution path)
mpif90=$MPIF90     # Fortran 90 parallel compiler (must be in execution path)
fflags=$FFLAGS     # Flags for Fortran 77 and 90 compilers
fflags_nomain=$FFLAGS_NOMAIN # Flags for linking Fortran sources with main in a different language
fflags_noopt=$FFLAGS_NOOPT # as FFLAGS With optimization disabled
f90flags=$F90FLAGS # Flags for Fortran 90 compiler only
ld=$LD             # Loader (must be in the execution path)
ldflags=$LDFLAGS   # Flags for loader
ld_libs=$LD_LIBS   # Additional libraries
blas_libs=$BLAS_LIBS     # blas library - specify e.g. /my/blas/lib/libmyblas.a
                         # or -L/my/blas/lib -lmyblas
lapack_libs=$LAPACK_LIBS # lapack library, similar to above
fft_libs=$FFT_LIBS       # FFT libraries - may depend upon DFLAGS
mpi_libs=$MPI_LIBS       # MPI libraries - shouldn't be needed
mass_libs=$MASS_LIBS     # MASS libraries (IBM only)
libdirs=$LIBDIRS         # Where to look for libraries (e.g. /my/blas/lib)
scalapack_libs=$SCALAPACK_LIBS # scalapack libs
scalapack_dir=$SCALAPACK_LIB  # Where to look for scalapack libs
blacs_dir=$BLACS_LIB          # Where to look for libblacs.a
ar=$AR                   # ar (shouldn't be needed)
arflags=$ARFLAGS         # Flags for ar (as above)
extlib_flags=$EXTLIB_FLAGS # Flags for internal copies of lapack and blas

  ]
)
