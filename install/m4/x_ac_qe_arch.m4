# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_ARCH], [

  AC_MSG_CHECKING([ARCH])

# many HPC systems are configured so that running parallel programs
# interactively is disabled: on those systems, AC_PROG_FC / _F77 / _CC
# would fail because they can't run the compiled executables.
# to work around that, let's pretend we are cross-compiling even if we aren't
# !!! this relies on undocumented Autoconf behavior !!!

# This is used to distinguish between true and fake cross compilation
# (only on NEC SX8 actually)
if test "$host" != "" ; then ranlib=echo; fi


# cross compiling? Why?
#cross_compiling=yes

if test "$host" = "" ; then host=$build; fi

# identify host architecture
if test "$arch" = ""
then
        case $host in
                ia64-*-linux-gnu )      arch=ia64   ;;
                x86_64-*-linux-gnu )    arch=x86_64 ;;
                arm-*linux* )           arch=arm    ;;
                *-pc-linux-gnu )        arch=ia32   ;;
                *-ibm-aix* )            arch=aix    ;;
                sparc-sun-* )           arch=sparc  ;;
                i386-pc-solaris* )      arch=solaris;;
                *86-apple-darwin* )     arch=mac686 ;;
                *-apple-darwin* )       arch=mac686 ;;
                *-pc-cygwin )           arch=cygwin ;;
                sx*-nec* )              arch=necsx  ;;
                powerpc64-*-linux-gnu ) arch=ppc64  ;;
                *-*-mingw32 )           arch=mingw32;;
                *-*-mingw64 )           arch=mingw64;;
                * )                     AC_MSG_WARN(Unrecognized build architecture)
        ;;
        esac
            # workaround for Cray-XT machines
        test -d /proc/cray_xt && arch=crayxt
            # workaround for IBM BG machines
        test -d /bgsys && arch=ppc64-bg
        test -f /bgsys/drivers/ppcfloor/bin/runjob && arch=ppc64-bgq
fi

  AC_MSG_RESULT(${arch})
  AC_SUBST(arch)

])
