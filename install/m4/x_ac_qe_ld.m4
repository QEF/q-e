# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_LD], [

# linker and archiver
# note that from this point on, further additions to
# linker flags should be added to ldflags rather than try_ldflags
if test "$ld" = ""       ; then ld="$mpif90"           ; fi
if test "$ldflags" = ""  ; then ldflags="$try_ldflags" ; fi
echo setting LD... $ld
echo setting LDFLAGS... $ldflags

# compilation flags for all subsequent tests
test_ldflags="`echo $ldflags | sed 's/\$([[^)]]*)//g'`"

AC_SUBST(ld)
AC_SUBST(ldflags)
  
  ]
)
