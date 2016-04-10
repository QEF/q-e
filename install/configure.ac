# Copyright (C) 2001-2016 Quantum ESPRESSO group
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

AC_INIT(ESPRESSO, 5.3.0, , espresso)

AC_PREREQ(2.60)
AC_CONFIG_MACRO_DIR([m4/])

# Initialize variables, filling with one comes from the environment...
X_AC_QE_DEFAULT_ENV()

# configure for current directory by default
if test "$topdir" = "" ; then topdir="`pwd`" ; fi

# check system type (no cross-compilation for now)
AC_CANONICAL_BUILD

# Checking Architecture...
X_AC_QE_ARCH()

# Add all needed -D options to try_dflags
try_dflags=""
try_iflags="-I../include"

# Checking archiver...
X_AC_QE_AR()

# Checking MPIF90...
X_AC_QE_MPIF90()

# Check environ...
X_AC_QE_ENVIRON()

# Checking MPIF90...
X_AC_QE_CC()

# Checking FFTW pointers...
X_AC_QE_FFTW_POINTER()

# Checking F77...
X_AC_QE_F77()

# Checking F90...
X_AC_QE_F90()

# Checking preprocessor...
X_AC_QE_CPP()

# Checking linker...
X_AC_QE_LD()

# Checking F90 rule...
X_AC_QE_F90RULE()

# find Fortran to C wrappers
X_AC_QE_WRAPPERS()

AC_LANG_PUSH(Fortran 77)
F77=$f90 # use Fortran 90 actually
FFLAGS="$test_fflags"
LDFLAGS="$test_ldflags"

# Checking BLAS...
X_AC_QE_BLAS()

# Checking LAPACK...
X_AC_QE_LAPACK()

# Checking for FFT...
X_AC_QE_FFT()

# Checking for IBM MASS library...
X_AC_QE_MASS()

# check for MPI library...
X_AC_QE_MPI()

# Setting ScaLAPACK...
X_AC_QE_SCALAPACK()

# Setting ELPA...
X_AC_QE_ELPA()

# Setting HDF5...
#X_AC_QE_HDF5()

# Checking SIGNAL...
X_AC_QE_SIGNAL()

# xlf compilers (AIX and powerpc) want comma-separated -D directives
X_AC_QE_AIX_DFLAGS()

# Checking for ranlib...
X_AC_QE_RANLIB()

# Check if wget or curl work  
X_AC_QE_WGET()

if test "$dflags" = "" ; then dflags="$try_dflags" ; fi
echo setting DFLAGS... $dflags

if test "$iflags" = "" ; then iflags="$try_iflags" ; fi
echo setting IFLAGS... $iflags

# export additional settings to generated files
AC_SUBST(dflags)
AC_SUBST(fdflags)
AC_SUBST(iflags)
AC_SUBST(ld_libs)
AC_SUBST(topdir)
AC_SUBST(extlib_flags)

AC_CONFIG_HEADERS([include/c_defs.h:include/c_defs.h.in])
AC_CONFIG_FILES([include/fft_defs.h:include/fft_defs.h.in])
AC_CONFIG_FILES([include/configure.h:include/configure.h.in])

#AC_CONFIG_FILES(Makefile.cpu)
AC_CONFIG_FILES(make.sys)
AC_CONFIG_FILES(configure.msg)
AC_CONFIG_FILES(install/make_wannier90.sys)


AC_OUTPUT

# final messages
sed '/@delete@/d' configure.msg
echo configure: success
