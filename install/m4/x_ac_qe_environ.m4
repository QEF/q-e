# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation
#
AC_DEFUN([X_AC_QE_ENVIRON], [

	AC_ARG_ENABLE(environment,
		[AS_HELP_STRING([--enable-environment], [compile solvent-related stuff (default: no)])],
		[if   test "$enableval" = "yes" ; then
			enable_environment=1
		else
			enable_environment=0
		fi],
		[enable_environment=0])
	
	if test "$enable_environment" -eq 1 ;
	then
		try_dflags="$try_dflags -D__ENVIRONMENT"
	fi  
	
	]
)
