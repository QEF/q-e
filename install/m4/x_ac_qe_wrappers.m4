# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_WRAPPERS], [

	AC_ARG_ENABLE(wrappers,
		[AS_HELP_STRING([--disable-wrappers], [disable C to fortran wrapper check (default: enabled)])],
		[if   test "$enableval" = "yes" ; then
			check_wrappers=1
		else
			check_wrappers=0
		fi],
		[check_wrappers=1])
   
	# Find Fortran to C wrappers
	if test "$check_wrappers" -ne 0; then
		AC_F77_WRAPPERS
	fi

	]
)
