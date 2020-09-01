# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_BEEF], [

beef_libs_switch=

AC_MSG_CHECKING([BEEF])

beef_libs=$BEEF_LIBS

AC_ARG_WITH(libbeef,
   [AS_HELP_STRING([--with-libbeef],
       [(yes|no) link against library for BEEF xc support (default: yes)])],
   [if   test "$withval" = "no" ; then
      with_beef=0
   else
      with_beef=1
   fi],
   [with_beef=1])

AC_ARG_WITH(libbeef-prefix,
   [AS_HELP_STRING([--with-libbeef-prefix],
       [path to directory containing libbeef.a (default: none)])],
   [libbeef_prefix="$withval"
   ],
   [libbeef_prefix=""])

if test "$with_beef" -ne 0
then
    BEEF_LDFLAGS_SAVE="$LDFLAGS"
    LDFLAGS="$beef_libs $LDFLAGS"
    if ! test "$libbeef_prefix"x = "x"
    then
        beefprefix=1
        LDFLAGS="$LDFLAGS -L$libbeef_prefix/lib -L$libbeef_prefix"

        LIBS="$LIBS -lbeef"
        AC_LINK_IFELSE([
          SUBROUTINE ddot
          END SUBROUTINE
          SUBROUTINE dgemv
          END SUBROUTINE
          PROGRAM main
          CALL beefx
          END
                       ], [usebeef=1], usebeef=0)
    else
        beefprefix=0
        usebeef=1
    fi

    if test "$usebeef" -ne 0
    then
        if test "$beefprefix" -ne 0
        then
            beef_libs="$beef_libs -L$libbeef_prefix/lib -L$libbeef_prefix -lbeef"
            beef_libs_switch="external"
        else
            beef_libs="\$(TOPDIR)/LIBBEEF/libbeef.a"
            beef_libs_switch="internal"
        fi
        echo "-lbeef"
        echo setting BEEF_LIBS... $beef_libs
        AC_SUBST(beef_libs)
    else
        echo no
        try_dflags="$try_dflags -D__NOBEEF"
    fi
    LDFLAGS="$BEEF_LDFLAGS_SAVE"
else
    echo no
    try_dflags="$try_dflags -D__NOBEEF"
fi

AC_SUBST(beef_libs_switch)
]

)
