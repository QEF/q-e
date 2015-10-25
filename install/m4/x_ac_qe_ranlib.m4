# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation
#
AC_DEFUN([X_AC_QE_RANLIB], [

 
  if test "$ranlib" != "echo"
  then
    AC_CHECK_PROG(ranlib,ranlib,ranlib,echo)
  fi

  if test "$arch" = "mac686"; then
    if test "$ranlib" = "ranlib"; then
      ranlib="ranlib -c"
    fi
  fi
  
  AC_SUBST(ranlib)
  
  ]
)
