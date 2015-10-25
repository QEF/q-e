# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_SIGNAL], [

   AC_ARG_ENABLE(signals,
     [AS_HELP_STRING([--enable-signals],
         [enable signal trapping (default: no)])],
     [if   test "$enableval" = "yes" ; then
        use_signals=1
     else
        use_signals=0
     fi],
     [use_signals=0])
   
  # preprocessing flag for signal trapping (experimental)
  if test "$use_signals" -eq 1 ; then try_dflags="$try_dflags -D__TRAP_SIGUSR1" ; fi
  
  ]
)
