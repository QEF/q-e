# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_BEEF], [

AC_MSG_CHECKING([BEEF])

beef_libs=$BEEF_LIBS

AC_ARG_ENABLE(libbeef,
   [AS_HELP_STRING([--enable-libbeef],
       [link against library for BEEF xc support (default: yes - if found)])],
   [if   test "$withval" = "yes" ; then
      beef_enabled=1
   else
      beef_enabled=0
   fi],
   [beef_enabled=1])

AC_ARG_WITH(libbeef-prefix,
   [AS_HELP_STRING([--with-libbeef-prefix],
       [path to directory containing libbeef.a (default: none)])],
   [libbeef_prefix="$withval"
   ],
   [libbeef_prefix=""])

if test "$beef_enabled" -ne 0
then
    BEEF_LDFLAGS_SAVE="$LDFLAGS"
    LDFLAGS="$beef_libs $LDFLAGS"
    if ! test "$libbeef_prefix"x = "x"
    then
        beefprefix=1
        LDFLAGS="$LDFLAGS -L$libbeef_prefix/lib -L$libbeef_prefix"
    else
        beefprefix=0
    fi
    LIBS="$LIBS -lbeef"
    AC_LINK_IFELSE([
      SUBROUTINE ddot
      END SUBROUTINE
      SUBROUTINE dgemv
      END SUBROUTINE
      PROGRAM main
      CALL beefx
      END
                   ], [try_dflags="$try_dflags -Duse_beef" ; usebeef=1], usebeef=0)
    if test "$usebeef" -ne 0
    then
        if test "$beefprefix" -ne 0
        then
            beef_libs="$beef_libs -L$libbeef_prefix/lib -L$libbeef_prefix -lbeef"
        else
            if test "$beef_libs"x = "x"
            then
                beef_libs=-lbeef
            fi
        fi
        echo "-lbeef"
        echo setting BEEF_LIBS...
        AC_SUBST(beef_libs)
    else
        echo no
    fi
    LDFLAGS="$BEEF_LDFLAGS_SAVE"
fi

]

)
