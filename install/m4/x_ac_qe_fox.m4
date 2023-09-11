# Copyright (C) 2022 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_FOX], [

   AC_ARG_WITH(fox,
     [AS_HELP_STRING([--with-fox],
         [use FoX library (default: no)])],
     [if   test "$withval" = "yes" ; then
        with_fox=1
     else
        with_fox=0
     fi],
     [with_fox=0])
   
  # Use internal code for xml read/write
  if test "$with_fox" -eq 0 ; then
     extfox=""
     foxmods=""
     foxlibs=""
  else
     extfox="libfox"
     foxmods="\$(MOD_FLAG)\$(TOPDIR)/FoX/finclude"
     foxlibs="-L\$(TOPDIR)/FoX/lib  -lFoX_dom -lFoX_sax -lFoX_wxml -lFoX_common -lFoX_utils -lFoX_fsys "
     try_iflags="$try_iflags -I\$(TOPDIR)/FoX/finclude "
     try_dflags="$try_dflags -D__fox"
fi

AC_SUBST(extfox)
AC_SUBST(foxmods)
AC_SUBST(foxlibs)

  ]
)
