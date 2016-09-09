# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_MPIF90], [

AC_ARG_ENABLE(parallel,
   [AS_HELP_STRING([--enable-parallel],
       [compile for parallel execution if possible (default: yes)])],
   [set_use_parallel=1
    if   test "$enableval" = "yes" ; then
      use_parallel=1
   else
      use_parallel=0
   fi],
   [set_use_parallel=0 use_parallel=1])
   
# candidate fortran compilers good for all cases
try_mpif90="mpif90"
try_f90="gfortran g95 f90"

# candidate compilers and flags based on architecture
case $arch in
ia32 | ia64 | x86_64 )
        try_f90="ifort pgf90 pathf95 sunf95 openf95 nagfor $try_f90"
        ;;
arm )
        try_f90="$try_f90"
        ;;
solaris )
        try_f90="sunf95 $try_f90"
        ;;
aix )
        try_mpif90="mpxlf90_r mpxlf90"
        try_f90="xlf90_r xlf90 $try_f90"
        try_dflags="-D__AIX -D__XLF"
        ;;
crayxt* )
        try_f90="ftn"
        try_mpif90="ftn"
        ;;
mac686 | cygwin )
        try_f90="ifort $try_f90"
        ;;
mingw* )
        ld="$F90"
        # this is set for C/C++, but we need it for Fortran, too.
        try_dflags="-D_WIN32"
        ;;
necsx )
        # most likely the following generates a bug
        sxopt=`echo $host|awk '{print substr($1,1,3)}'`
        echo $sxopt $host
        try_mpif90="sxmpif90"
        try_f90="sxf90"
        try_dflags='-D__SX6 '
        use_fft_asl=0
        use_fft_mathkeisan=1
        use_fft_para=0
# default for Nec: no parallel unless explicitly required
        if test "$set_use_parallel" -ne 1 ; then use_parallel=0 ; fi
        if test "$use_parallel" -eq 1 ; then use_fft_para=1 ; fi
        try_dflags_fft_asl='-DASL'
        try_dflags_fft_mathkeisan=' '
        try_dflags_fft_para='-D__USE_3D_FFT'
        ;;
ppc64 )
        try_mpif90="mpxlf90_r mpf90_r mpif90"
        try_f90="xlf90_r $try_f90"
        try_dflags="-D__XLF"
        ;;
# PowerPC MareNostrum
ppc64-mn )
        try_f90="xlf90_r"
        try_dflags="-D__XLF"
        ;;
# IBM BlueGene
ppc64-bg | ppc64-bgq )
	if test "$use_openmp" -eq 0 ; then
          try_mpif90="mpixlf90"
          try_f90="bgxlf90"
	else
          try_mpif90="mpixlf90_r"
          # Executable paths are usually consistent across several 
          # IBM BG/P BG/Q machine deployed 
          ld="/bgsys/drivers/ppcfloor/comm/xl.ndebug/bin/mpixlf90_r"
          try_f90="bgxlf90_r"
	fi
        try_arflags="ruv"
        try_dflags="-D__XLF"
        ;;
* )
        AC_MSG_WARN($arch : unsupported architecture?)
        ;;
esac

# check serial Fortran 90 compiler. This must be done before performing
# the check for the parallel compiler (section below) because option
# --disable-parallel will do very strange things otherwise. The reason
# seems to be that autoconf does not repeat all tests for the second
# occurrence of AC_PROG_FC. So the first occurrence is the one that
# must always be performed, the second is optional. PG & CC sep.2006

# use F90 if set
if test "$f90" = "" ; then f90="$try_f90" ; fi
AC_PROG_FC($f90)
f90=$FC
AC_FC_SRCEXT(f90)

# check parallel Fortran 90 compiler
if test "$use_parallel" -eq 0 ;
then
        mpif90=$f90
else
        # clear cached values (not sure when and why this is needed)
        unset FC ac_cv_prog_ac_ct_FC ac_cv_fc_compiler_gnu ac_cv_prog_fc_g
        if test "$mpif90" = "" ; then 
	   mpif90="$try_mpif90 $f90"
           AC_PROG_FC($mpif90)
        else
           AC_PROG_FC($mpif90)
# this avoids that an empty MPIF90 field is produced if the corresponding
# environment variable MPIF90 does not contain an acceptable compiler
           if test "$FC" = "" ; then 
		AC_MSG_WARN([MPIF90 not found: using MPIF90 anyway])
	  	FC=$MPIF90
	   fi
        fi
        mpif90=$FC
fi

# check which compiler does mpif90 wrap

case "$arch" in
        ia32 | ia64 | x86_64 | mac686 )
        echo $ECHO_N "checking version of $mpif90... $ECHO_C"
        ifort_version=`$mpif90 -V 2>&1 | grep "Intel(R)"`
        sunf95_version=`$mpif90 -V 2>&1 | grep "Sun Fortran"`
        openf95_version=`$mpif90 -V 2>&1 | grep "^Open64"`
        pgf_version=`$mpif90 -V 2>&1 | grep "^pgf"`
        g95_version=`$mpif90 -v 2>&1 | grep "g95"`
        enzo_version=`$mpif90 -v 2>&1 | grep "PathScale ENZO"`
        eko_version=`$mpif90 -v 2>&1 | grep "PathScale EKOPath"`
        pathf95_version=`$mpif90 -v 2>&1 | grep "PathScale"`
        gfortran_version=`$mpif90 -v 2>&1 | grep "gcc version"`
        nagfor_version=`$mpif90 -v 2>&1 | grep "NAG Fortran"`
        #
        if test "$ifort_version" != ""
        then
                version=`$mpif90 --version 2>&1 | grep "IFORT" | cut -d ' ' -f3`
                f90_major_version=`echo $version | cut -d. -f1`
                echo "${ECHO_T}ifort $f90_major_version"
                f90_in_mpif90="ifort"
                # Why so?
                if test "$f90_major_version" -gt "9"; then
                   MKL_FLAGS="-static-intel"
                fi

        elif test "$sunf95_version" != ""
        then
                version=`echo $sunf95_version | cut -d ' ' -f5`
                echo "${ECHO_T}sunf95 $version"
                f90_in_mpif90="sunf95"
        elif test "$openf95_version" != ""
        then
                version=`echo $openf95_version | cut -d ' ' -f5`
                echo "${ECHO_T}openf95 $version"
                f90_in_mpif90="openf95"
        elif test "$pgf_version" != ""
        then
                version=`echo $pgf_version | cut -d ' ' -f2`
                echo "${ECHO_T}pgf90 $version"
                f90_in_mpif90="pgf90"
                # flag to test MKL with PGI
                MKL_FLAGS="-pgf90libs"
        elif test "$enzo_version" != ""
        then
                version=`echo $enzo_version | cut -d ' ' -f6`
                echo "${ECHO_T}pathf95 $version"
                f90_in_mpif90="pathf95"
        elif test "$eko_version" != ""
        then
                version=`echo $eko_version | cut -d ' ' -f6`
                echo "${ECHO_T}pathf95 $version"
                f90_in_mpif90="pathf95"
        elif test "$g95_version" != ""
        then
                version=`echo $g95_version | cut -d ' ' -f3`
                echo "${ECHO_T}g95 $version"
                f90_in_mpif90="g95"
        elif test "$pathf95_version" != ""
        then
                version=`echo $pathf95_version | cut -d ' ' -f5`
                echo "${ECHO_T}pathf95 $version"
                f90_in_mpif90="pathf95"
        elif test "$gfortran_version" != ""
        then
                version=`echo $gfortran_version | cut -d ' ' -f3`
                f90_major_version=`echo $version | cut -d. -f1`
                f90_minor_version=`echo $version | cut -d. -f2`
                echo "${ECHO_T}gfortran $f90_major_version.$f90_minor_version"
                f90_in_mpif90="gfortran"
        elif test "$nagfor_version" != ""
        then
                # NAG 6.0 has the codename attached to version number... annoying
                version=`echo $nagfor_version | cut -d ' ' -f5`
                echo "${ECHO_T}nagfor $version"
                f90_in_mpif90="nagfor"
        else
                echo "${ECHO_T}unknown, assuming gfortran"
                f90_in_mpif90="gfortran"
        fi
        # check if serial and parallel compiler are the same
        if test "$f90" != "$f90_in_mpif90"; then
           AC_MSG_WARN([parallel compiler $mpif90 uses $f90_in_mpif90, but serial compiler $f90 was detected])
           AC_MSG_WARN([assuming F90=$f90_in_mpif90, discarding $f90])
        fi
        f90=$f90_in_mpif90
        ;;
esac

echo setting F90... $f90
echo setting MPIF90... $mpif90

case "$f90" in
f90 | fc | ftn )
    echo $ECHO_N "checking version wrapped by $f90 command... $ECHO_C"

    if $f90 -V 2>&1 | grep -q "Intel(R)" ; then
        f90_flavor=ifort
    elif $f90 -V 2>&1 | grep -q "Sun Fortran" ; then
        f90_flavor=sunf95
    elif $f90 -V 2>&1 | grep -q "^Open64" ; then
        f90_flavor=openf95
    elif $f90 -V 2>&1 | grep -q "^pgf" ; then
        f90_flavor=pgf
    elif $f90 -v 2>&1 | grep -q "PathScale ENZO" ; then
        f90_flavor=pathf95
    elif $f90 -v 2>&1 | grep -q "PathScale EKOPath" ; then
        f90_flavor=pathf95
    elif $f90 -version 2>&1 | grep -q "PathScale" ; then
        f90_flavor=pathf95
    elif $f90 -v 2>&1 | grep -q "g95" ; then
        f90_flavor=g95
    elif $f90 -v 2>&1 | grep -q "gcc version" ; then
        f90_flavor=gfortran
    elif $f90 -V 2>&1 | grep -q "Cray Fortran" ; then
        f90_flavor=crayftn
    elif $f90 -version 2>&1 | grep -q "NAG Fortran" ; then
        f90_flavor=nagfor
    else
        echo $ECHO_N "unknown, leaving as... $ECHO_C"
        f90_flavor=$f90
    fi
    echo $f90_flavor
    ;;
* )
    f90_flavor=$f90
    ;;
esac

AC_SUBST(f90)
AC_SUBST(mpif90)

])
