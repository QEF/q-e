# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_ELPA], [

  AC_MSG_CHECKING([ELPA])
 
AC_ARG_WITH(elpa-include,
   [AS_HELP_STRING([--with-elpa-include],
       [Specify full path ELPA include and modules headers (default: no)])],
   [if  test "$withval" = "no" ; then
      with_elpa_include=0
   else  
      with_elpa_include=1
      elpa_include="$withval"
   fi],
   [with_elpa_include=0])

AC_ARG_WITH(elpa-lib,
   [AS_HELP_STRING([--with-elpa-lib],
       [Specify full path ELPA static or dynamic library (default: no)])],
   [if  test "$withval" = "no" ; then
      with_elpa_libs=0
   else
      with_elpa_libs=1
      elpa_libs="$withval"
   fi],
   [with_elpa_libs=0])


elpa_line="@delete@"

# ELPA iff SCALAPACK
if test "$with_elpa_libs" -eq 1; then 
  if test "$have_scalapack" -eq 1; then
    if test "$with_elpa_include" -eq 1 && test "$with_elpa_libs" -eq 1; then
      scalapack_libs="$elpa_libs $scalapack_libs"
      try_iflags="$try_iflags -I$elpa_include"
      try_dflags="$try_dflags -D__ELPA_2016"
      elpa_line="ELPA_LIBS=$elpa_libs"
    fi
  else
    AC_MSG_WARN([*** ScaLAPACK is needed to use ELPA])
  fi
fi

  AC_MSG_RESULT(${elpa_libs})
 
  AC_SUBST(elpa_libs)
  AC_SUBST(elpa_line) 
  ]
)
