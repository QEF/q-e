# Copyright (C) 2001-2025 Quantum ESPRESSO Foundation

# Files (mostly Makefiles) to be copied to the build directories
# NB: make.depend files and Makefiles including them are copied in makedeps.sh, not here

AC_DEFUN([X_AC_QE_MAKE_TREE], [

AC_CONFIG_FILES([install/extlibs_makefile:install/extlibs_makefile])
AC_CONFIG_FILES([install/install_utils:install/install_utils])
AC_CONFIG_FILES([install/plugins_makefile:install/plugins_makefile])
AC_CONFIG_FILES([install/plugins_list:install/plugins_list])
AC_CONFIG_FILES([install/tldeps:install/tldeps])

AC_CONFIG_FILES([Makefile:Makefile])

AC_CONFIG_FILES([COUPLE/Makefile:COUPLE/Makefile])
AC_CONFIG_FILES([COUPLE/src/Makefile:COUPLE/src/Makefile])

AC_CONFIG_FILES([FFTXlib/Makefile:FFTXlib/Makefile])
AC_CONFIG_FILES([FFTXlib/tests/Makefile:FFTXlib/tests/Makefile])

AC_CONFIG_FILES([PIOUD/Makefile:PIOUD/Makefile])
AC_CONFIG_FILES([PIOUD/src/Makefile:PIOUD/src/Makefile])
AC_CONFIG_FILES([PIOUD/src/make.depend:PIOUD/src/make.depend])

AC_CONFIG_FILES([NEB/Makefile:NEB/Makefile])

AC_CONFIG_FILES([XSpectra/Makefile:XSpectra/Makefile])

AC_CONFIG_FILES([XClib/Makefile:XClib/Makefile])

AC_CONFIG_FILES([Modules/Makefile:Modules/Makefile])

AC_CONFIG_FILES([atomic/Makefile:atomic/Makefile])

AC_CONFIG_FILES([PP/Makefile:PP/Makefile])
AC_CONFIG_FILES([PP/simple_transport/src/Makefile:PP/simple_transport/src/Makefile])

AC_CONFIG_FILES([QEHeat/Makefile:QEHeat/Makefile])

AC_CONFIG_FILES([GWW/Makefile:GWW/Makefile])
AC_CONFIG_FILES([GWW/util/Makefile:GWW/util/Makefile])
AC_CONFIG_FILES([GWW/minpack/Makefile:GWW/minpack/Makefile])

AC_CONFIG_FILES([HP/Makefile:HP/Makefile])

AC_CONFIG_FILES([PW/Makefile:PW/Makefile])

AC_CONFIG_FILES([PWCOND/Makefile:PWCOND/Makefile])

AC_CONFIG_FILES([KS_Solvers/Makefile:KS_Solvers/Makefile])

AC_CONFIG_FILES([CPV/Makefile:CPV/Makefile])

# TG: TDDFPT/ColorCalculator/Install?
AC_CONFIG_FILES([TDDFPT/Makefile:TDDFPT/Makefile])

# TG: autotest.inc?
AC_CONFIG_FILES([UtilXlib/tests/Makefile:UtilXlib/tests/Makefile])

AC_CONFIG_FILES([LAXlib/tests/Makefile:LAXlib/tests/Makefile])

AC_CONFIG_FILES([PHonon/Makefile:PHonon/Makefile])

AC_CONFIG_FILES([KCW/Makefile:KCW/Makefile])

AC_CONFIG_FILES([EPW/Makefile:EPW/Makefile])
AC_CONFIG_FILES([EPW/irobjs/Makefile:EPW/irobjs/Makefile])

AC_CONFIG_FILES([archive/README.md:archive/README.md])

]
)
