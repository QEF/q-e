# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_AIX_DFLAGS], [

# xlf compilers (AIX and powerpc) want comma-separated -D directives
if test "$xlf_flags" -ne 0
then
        fdflags="`echo $dflags | sed 's/  */,/g'`"
else
        fdflags="\$(DFLAGS)"
fi
  
  ]
)
