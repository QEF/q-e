-# Copyright (C) 2001-2023 Quantum ESPRESSO Foundation
 
AC_DEFUN([X_AC_QE_ENABLE_OPTIONS], [

AC_ARG_ENABLE(legacy_plugins,
	[ --enable-legacy_plugins    Enables use of legacy plugins in PW and CP.],
	[enable_legacy_plugins=$enableval],
	[enable_legacy_plugins="no"])
AC_MSG_RESULT(checking if legacy plugins are enabled ...  $enable_legacy_plugins) 
if test "$enable_legacy_plugins" = "yes";then
	try_dflags="$try_dflags -D__LEGACY_PLUGINS"
fi 
]
)
