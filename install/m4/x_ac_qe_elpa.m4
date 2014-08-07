# Copyright (C) 2001-2014 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_ELPA], [

  AC_MSG_CHECKING([ELPA])
 
AC_ARG_WITH(elpa,
   [AS_HELP_STRING([--with-elpa],
       [(yes|no|<path>) Use ELPA. Self-compile or a <path> can be specified (default: no)])],
   [if  test "$withval" = "yes" ; then
      with_elpa=1
   elif  test "$withval" = "no" ; then
      with_elpa=0
   else
      with_elpa=2
      with_elpa_path="$withval"
   fi],
   [with_elpa=0])

elpa_libs=""
 
# ELPA iff SCALAPACK (statically linked)
elpa_libs_switch="disabled"
if test "$with_elpa" -eq 2 && test "$have_scalapack" -eq 1; then
    if test "$use_openmp" -eq 1 ; then
        elpa_libs="$with_elpa_path/lib/libelpa_mt.a"
        try_iflags="$try_iflags -I$with_elpa_path/include/elpa/modules "
        try_dflags="$try_dflags -D__ELPA"
    else
        elpa_libs="$with_elpa_path/lib/libelpa.a"
        try_iflags="$try_iflags -I$with_elpa_path/include/elpa/modules "
        try_dflags="$try_dflags -D__ELPA"
    fi
  scalapack_libs="$elpa_libs $scalapack_libs"
fi

if test "$with_elpa" -eq 1 && test "$have_scalapack" -eq 1; then
    elpa_libs="\$(TOPDIR)/ELPA/libelpa.a"
    scalapack_libs="$elpa_libs $scalapack_libs"
    try_dflags="$try_dflags -D__ELPA"
    elpa_libs_switch="enabled"
fi

  AC_MSG_RESULT(${elpa_libs})
  AC_SUBST(elpa_libs_switch)
  
  ]
)
