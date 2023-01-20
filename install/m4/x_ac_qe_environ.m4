# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation
#
AC_DEFUN([X_AC_QE_ENVIRON], [

	AC_MSG_CHECKING([for Environ])

	AC_ARG_WITH(environ,
		[AS_HELP_STRING([--with-environ], [absolute path to Environ root directory])],
		[
			if test "$withval" != yes ;
		 	then
		 	 	environ_root="$withval"
		 	else
				environ_root="$topdir/Environ"
		 	fi
		 	
			if test -d $environ_root; then
				environ_libs="-L$environ_root/libs -lenvsrc -lenvfft -lenvutil"
				try_iflags="$try_iflags -I$environ_root/src"
		 		try_dflags="$try_dflags -D__ENVIRON"
				AC_MSG_RESULT(found at $environ_root)
			else
			    AC_MSG_ERROR([$environ_root is not a valid path])
			fi

		], 
		[AC_MSG_RESULT(not used)]
	)

	AC_SUBST(environ_libs)
	]
)
