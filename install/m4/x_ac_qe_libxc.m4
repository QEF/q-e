## Copyright (C) 2010-2015 M. Marques, X. Andrade, D. Strubbe, M. Oliveira
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##
## $Id: libxc.m4 12311 2016-04-19 20:52:52Z dstrubbe $
##

## Modified by D. Ceresoli to work with Quantum-Espresso

AC_DEFUN([ACX_LIBXC], [
acx_libxc_ok=no

AC_ARG_WITH(libxc,
   [AS_HELP_STRING([--with-libxc],
       [(yes|no) Use libXC for some XC functionals (default: no)])],
   [if  test "$withval" = "no" ; then
      with_libxc=0
   else
      with_libxc=1
   fi],
   [with_libxc=0]
)

AC_ARG_WITH(libxc-prefix, [AS_HELP_STRING([--with-libxc-prefix=DIR], [Directory where libxc was installed.])])
AC_ARG_WITH(libxc-include, [AS_HELP_STRING([--with-libxc-include=DIR], [Directory where libxc Fortran headers were installed.])])

dnl continue only if --with-libxc=yes
if test "$with_libxc" -ne 0; then

dnl Set FCFLAGS_LIBXC only if not set from environment
if test x"$FCFLAGS_LIBXC" = x; then
  case $with_libxc_prefix in
    "") FCFLAGS_LIBXC="$imod/usr/include" ;;
    *)  FCFLAGS_LIBXC="$imod$with_libxc_prefix/include" ;;
  esac
fi

case $with_libxc_include in
  "") ;;
  *)  FCFLAGS_LIBXC="$imod$with_libxc_include" ;;
esac

dnl Backup LIBS and FCFLAGS
acx_libxc_save_LIBS="$LIBS"
acx_libxc_save_FCFLAGS="$FCFLAGS"

dnl The tests
AC_MSG_CHECKING([for libxc])
AC_LANG_PUSH(Fortran)
FCFLAGS="$FCFLAGS_LIBXC $acx_libxc_save_FCFLAGS"

testprog="AC_LANG_PROGRAM([],[
  use xc_f90_lib_m
  implicit none
  integer :: major
  integer :: minor
  integer :: micro
  call xc_f90_version(major, minor, micro)])"

dnl set from environment variable, if not blank
if test ! -z "$LIBS_LIBXC"; then
  LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
  AC_LINK_IFELSE($testprog, [acx_libxc_ok=yes], [])
fi

if test ! -z "$with_libxc_prefix"; then
  if test x"$acx_libxc_ok" = xno; then
    LIBS_LIBXC="-L$with_libxc_prefix/lib -lxcf90 -lxc"
    LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
    AC_LINK_IFELSE($testprog, [acx_libxc_ok=yes], [])
  fi
else
  LIBS_LIBXC="-lxcf90 -lxc"
  LIBS="$LIBS_LIBXC $acx_libxc_save_LIBS"
  AC_LINK_IFELSE($testprog, [acx_libxc_ok=yes], [])
fi

AC_MSG_RESULT([$acx_libxc_ok ($FCFLAGS_LIBXC $LIBS_LIBXC)])

dnl Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
libxc_line="@delete@"
if test x"$acx_libxc_ok" = xyes; then
  try_dflags="$try_dflags -D__LIBXC"
  try_iflags="$try_iflags $FCFLAGS_LIBXC"
  AC_SUBST(LIBS_LIBXC)
  libxc_line="LIBXC_LIBS= $LIBS_LIBXC"
  AC_SUBST(libxc_line)
else
  AC_MSG_ERROR([Could not find required libxc library.])
fi

AC_LANG_POP(Fortran)

fi  dnl with_libxc=yes
])


