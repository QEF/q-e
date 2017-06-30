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
   
AC_ARG_WITH(elpa-version,
   [AS_HELP_STRING([--with-elpa-version],
       [Specify ELPA version, only year (default: 2016)])],
   [if  test "$withval" = "no" ; then
      with_elpa_version=0
   else
      with_elpa_version="$withval"
   fi],
   [with_elpa_version="2016"])


elpa_line="@delete@"

# ELPA iff SCALAPACK
if test "$with_elpa_libs" -eq 1; then 
  if test "$have_scalapack" -eq 1; then
    if test "$with_elpa_include" -eq 1 && test "$with_elpa_libs" -eq 1; then
      
      if test "$with_elpa_version" = "2015"; then
        try_dflags="$try_dflags -D__ELPA_2015"
      elif test "$with_elpa_version" = "2016"; then
        try_dflags="$try_dflags -D__ELPA_2016"
      elif test "$with_elpa_version" = "2017"; then
        try_dflags="$try_dflags -D__ELPA_2017"
      else
        AC_MSG_WARN([*** Invalid ELPA version, defaulting to 2016])
        try_dflags="$try_dflags -D__ELPA_2016"
      fi
      
      try_iflags="$try_iflags -I$elpa_include"
      scalapack_libs="$elpa_libs $scalapack_libs"
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
