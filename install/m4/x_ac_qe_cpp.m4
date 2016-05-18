# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_CPP], [

# preprocessor - try cpp in all cases; the preprocessor returned by
# AC_PROG_CPP
# may sometimes refuse to preprocess fortran files
if test "$cpp" = "" ; then cpp=$try_cpp; fi
# if test "$cpp" = "" ; then cpp=$CPP;     fi
echo setting CPP... $cpp

echo $ECHO_N "setting CPPFLAGS... $ECHO_C"
# Note: option -C makes trouble with recent gcc versions and pgi
case $cpp in
        cpp)  try_cppflags="-P -traditional" ;;
        fpp)  try_cppflags="-P "              ;;
        *)    try_cppflags=""                ;;
esac
if test "$cppflags" = "" ; then cppflags=$try_cppflags ; fi
echo "${ECHO_T}$cppflags"

# compilation flags for all subsequent tests
test_cppflags="$test_cflags"

AC_SUBST(cpp)
AC_SUBST(cppflags)
  
  ]
)
