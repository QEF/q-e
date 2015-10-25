# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_FFTW_POINTER], [

# check size of pointers to int - needed to decide the size of integer
# arrays in fortran holding C pointers for FFTW

AC_CHECK_SIZEOF([int *])
SIZEOF_INT_P=$ac_cv_sizeof_int_p
AC_SUBST(SIZEOF_INT_P)

# check if the structure mallinfo is present in malloc.h
AC_CHECK_HEADER(malloc.h,have_malloc_h=1,have_malloc_h=0, )
if test "$have_malloc_h" -ne 0
then
AC_CHECK_MEMBER([struct mallinfo.arena],
                [AC_DEFINE(HAVE_MALLINFO)],
                ,
                [#include <malloc.h>])

fi

])
