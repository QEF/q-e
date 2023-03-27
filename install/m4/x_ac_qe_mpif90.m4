# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_MPIF90], [
AC_REQUIRE([AC_PROG_FC])
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
try_f90="gfortran f90"

# candidate compilers and flags based on architecture
case $arch in
ia32 | ia64 | x86_64 )
        try_f90="ifort nvfortran pgf90 nagfor $try_f90"
        try_mpif90="mpiifort $try_mpif90"
        ;;
arm )
        try_f90="nvfortran pgf90 armflang $try_f90"
        ;;
craype* )
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
        ;;
# PowerPC little endian
ppc64le )
        try_mpif90="$try_mpif90 mpixlf"
        try_f90="xlf90_r"
        ;;
# IBM BlueGene - obsolete
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
        ;;
* )
        AC_MSG_WARN($arch : unsupported architecture?)
        ;;
esac

# check Fortran 90 compiler

# clear cached values
unset FC ac_cv_prog_ac_ct_FC ac_cv_fc_compiler_gnu ac_cv_prog_fc_g

if test "$use_parallel" -eq 0 ; then
# scalar case - use F90 if set
    	if test "$f90" = "" ; then
	   mpif90="$try_f90"
	else
	   mpif90="$f90"
	fi
else
# parallel case - use MPIF90 if set
        if test "$mpif90" = "" ; then 
	   mpif90="$try_mpif90 $f90 $try_f90 "
	fi
    	if test "$f90" != "" ; then
           AC_MSG_WARN([F90 value is set to be consistent with value of MPIF90])
	fi
fi

AC_PROG_FC($mpif90)
# this avoids that an empty MPIF90 field is produced if the corresponding
# environment variable MPIF90 does not contain an acceptable compiler
if test "$FC" = "" ; then 
   AC_MSG_WARN([MPIF90 not found: using MPIF90 anyway])
   FC=$mpif90
fi
mpif90=$FC

# check which compiler does mpif90 wrap

echo $ECHO_N "checking version of $mpif90... $ECHO_C"
ifort_version=`$mpif90 -V 2>&1 | grep "Intel(R)"`
pgf_version=`$mpif90 -V 2>&1 | grep "^pgf"`
nvfortran_version=`$mpif90 -V 2>&1 | grep "^nvfortran"`
gfortran_version=`$mpif90 -v 2>&1 | grep "gcc version"`
nagfor_version=`$mpif90 -v 2>&1 | grep "NAG Fortran"`
xlf_version=`$mpif90 -v 2>&1 | grep "xlf"`
armflang_version=`$mpif90 -v 2>&1 | grep "Arm C/C++/Fortran Compiler version"`
frt_version=`$mpif90 -v 2>&1 | grep "Fujitsu Fortran Compiler"`
ftn_version=`$mpif90 --version 2>&1 | grep "Cray Fortran"`
#
if test "$ifort_version" != ""
then
        version=`$mpif90 --version 2>&1 | grep "IFORT" | cut -d ' ' -f3`
        f90_major_version=`echo $version | cut -d. -f1`
        echo "${ECHO_T}ifort $f90_major_version"
        f90_in_mpif90="ifort"
elif test "$nvfortran_version" != ""
then
        version=`echo $nvfortran_version | cut -d ' ' -f2`
        f90_major_version=`echo $version | cut -d. -f1`
	f90_minor_version=`echo $version | cut -d. -f2 | cut -d- -f1`
        echo "${ECHO_T}nvfortran $version"
        f90_in_mpif90="nvfortran"
elif test "$pgf_version" != ""
then
        version=`echo $pgf_version | cut -d ' ' -f2`
        f90_major_version=`echo $version | cut -d. -f1`
        echo "${ECHO_T}pgf90 $version"
        f90_in_mpif90="pgf90"
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
elif test "$xlf_version" != ""
then
        echo "${ECHO_T}xlf (version unknonw)"
        f90_in_mpif90="xlf90_r"
        try_dflags="-D__XLF"
elif test "$armflang_version" != "" 
then 
        version=`echo $armflang_version | cut -d" " -f 5`
        f90_major_version=`echo $version | cut -d. -f1` 
        f90_minor_version=`echo $version | cut -d. -f2` 
        f90_in_mpif90="armflang"
        try_foxflags="-D__PGI"  
elif test "$frt_version" != ""
then
        version=`echo $frt_version | cut -d" " -f 5`
        echo "${ECHO_T}frt $version"
        f90_in_mpif90="frt"
        try_foxflags="-D__FUJITSU"
elif test "$ftn_version" != ""
then
        version=`echo $fnt_version | cut -d" " -f 5`
        echo "${ECHO_T}ftn $version"
        f90_in_mpif90="ftn"
else
        echo "${ECHO_T}unknown, assuming gfortran"
        f90_in_mpif90="gfortran"
fi
# notify if serial and parallel compiler are the same
if test "$set_use_parallel" -eq 1 ; then
   if test "$mpif90" = "$f90_in_mpif90" && test "$mpif90" != "ftn"; then
      AC_MSG_WARN([parallel and serial compiler are the same])
   fi
fi

case "$mpif90" in
ftn )
    # For cray compiler, mpif90=f90=ftn, f90_flavor used to select flags
    echo $ECHO_N "checking version wrapped by $mpif90 command... $ECHO_C"
    if $mpif90 -V 2>&1 | grep -q "Intel(R)" ; then
        f90_flavor=ifort
    elif $mpif90 -V 2>&1 | grep -q "^pgf" ; then
        f90_flavor=pgf
    elif $mpif90 -v 2>&1 | grep -q "gcc version" ; then
        f90_flavor=gfortran
    elif $mpif90 -V 2>&1 | grep -q "Cray Fortran" ; then
        f90_flavor=crayftn
    elif $mpif90 -version 2>&1 | grep -q "NAG Fortran" ; then
        f90_flavor=nagfor
    else
        echo $ECHO_N "unknown, leaving as ... $ECHO_C"
        f90_flavor=$mpif90
    fi
    echo $f90_flavor
    f90=ftn
    ;;
* )
    # For all other cases f90=f90_flavor=f90_in_mpif90
    f90=$f90_in_mpif90
    f90_flavor=$f90
    ;;
esac

echo setting F90... $f90
echo setting MPIF90... $mpif90

AC_SUBST(f90)
AC_SUBST(mpif90)

])
