# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_F77], [

# candidate f77 compilers good for all cases
try_f77="$f90"

case "$arch:$f90_flavor" in
*:pgf90 )
        try_f77="pgf77 $f90"
        ;;
aix:*xlf*_r )
        try_f77="xlf_r $f90"
        ;;
aix:*xlf* )
        try_f77="xlf $f90"
        ;;
sparc:* | solaris:sunf95 )
        try_f77="f77 $f90"
        ;;
ppc64-bg*:*xlf90_r )
        try_f77="bgxlf_r"
        ;;
ppc64-bg*:*xlf90 )
        try_f77="bgxlf"
        ;;
ppc64:*xlf* | ppc64-mn:*xlf* )
        try_f77="xlf_r $f90"
        ;;
esac

# check serial Fortran 77 compiler (use F77 if it was set)
if test "$f77" = "" ; then f77="$try_f77" ; fi
AC_PROG_F77($f77)
f77=$F77

echo setting F77... $f77

AC_SUBST(f77)

])
