AC_INIT(ESPRESSO, 5.1.svn, , espresso)
# ----------------------------
# Non-standard precious variables
# ----------------------------
AC_ARG_VAR(EXTLIB_FLAGS, This variable controls the flags passed to internal BLAS and LAPACK libraries)
# -----------------------------
# options
# -----------------------------

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

AC_ARG_ENABLE(openmp,
   [AS_HELP_STRING([--enable-openmp],
       [compile for openmp execution if possible (default: no)])],
   [if   test "$enableval" = "yes" ; then
      use_openmp=1
   else
      use_openmp=0
   fi],
   [use_openmp=0])

# debug flags are implemented only for a few cases
AC_ARG_ENABLE(debug,
   [AS_HELP_STRING([--enable-debug],
       [compile Fortran with debug flags (default: no)])],
   [if   test "$enableval" = "yes" ; then
      use_debug=1
   else
      use_debug=0
   fi],
   [use_debug=0])

# shared library flags are implemented only for a few (untested) cases
AC_ARG_ENABLE(shared,
   [AS_HELP_STRING([--enable-shared],
       [use shared libraries if available (default: yes)])],
   [if   test "$enableval" = "yes" ; then
      use_shared=1
   else
      use_shared=0
   fi],
   [use_shared=1])

# the  C-to-Fortran wrapper check may cause trouble on some machines
AC_ARG_ENABLE(wrappers,
   [AS_HELP_STRING([--disable-wrappers],
       [disable C to fortran wrapper check (default: enabled)])],
   [if   test "$enableval" = "yes" ; then
      check_wrappers=1
   else
      check_wrappers=0
   fi],
   [check_wrappers=1])

AC_ARG_ENABLE(signals,
   [AS_HELP_STRING([--enable-signals],
       [enable signal trapping (default: no)])],
   [if   test "$enableval" = "yes" ; then
      use_signals=1
   else
      use_signals=0
   fi],
   [use_signals=0])

# -----------------------------
# packages
# ----------------------------

AC_ARG_WITH(scalapack,
   [AS_HELP_STRING([--with-scalapack],
       [(yes|no|intel) Use scalapack if available. Set to "intel" to use Intel MPI and blacs (default: use openMPI)])],
   [if  test "$withval" = "yes" ; then
      with_scalapack=1
   elif  test "$withval" = "intel" ; then
      with_scalapack=2
   elif  test "$withval" = "no" ; then
      with_scalapack=0
   fi],
   [with_scalapack=1])

AC_ARG_WITH(elpa,
   [AS_HELP_STRING([--with-elpa],
       [(yes|no|<path>) Use ELPA. Self-compile or a <path> can be specified (default: no)])],
   [if  test "$withval" = "yes" ; then
      with_elpa=1
   elif  test "$withval" = "no" ; then
      with_elpa=0
   else
      with_elpa=2
      with_elpa_path="$withval"   
   fi],
   [with_elpa=0])

AC_ARG_WITH(internal-blas,
   [AS_HELP_STRING([--with-internal-blas],
       [compile with internal blas (default: no)])],
    [if   test "$withval" = "yes" ; then
      use_internal_blas=1
   else
      use_internal_blas=0
   fi],
   [use_internal_blas=0])

AC_ARG_WITH(internal-lapack,
   [AS_HELP_STRING([--with-internal-lapack],
       [compile with internal lapack (default: no)])],
    [if   test "$withval" = "yes" ; then
      use_internal_lapack=1
   else
      use_internal_lapack=0
   fi],
   [use_internal_lapack=0])
# -----------------------------
# more options 
# ----------------------------

AC_ARG_ENABLE(environment,
   [AS_HELP_STRING([--enable-environment],
       [compile solvent-related stuff (default: no)])],
    [if   test "$enableval" = "yes" ; then
      enable_environment=1
   else
      enable_environment=0
   fi],
   [enable_environment=0])

# store variables from the environment, if set (may be or not be set)
# If set, they take precedence over configure internal choice.
# Flags and libraries are accepted without further testing;
# compilers are tested. Specify compiler name only, not the full path
# (i.e. F90=/usr/local/bin/f90 may not work, use F90=f90)

topdir=$TOPDIR     # current directory
arch=$ARCH         # see below for recognized architectures
env_cc=$CC         # C compiler (must be in the execution path)
cpp=$CPP           # C preprocessor (as above)
cflags=$CFLAGS     # Flags for C compiler
cppflags=$CPPFLAGS # Flags for C preprocessor
dflags=$DFLAGS     # Fortran file preprocessing options, e.g. -D__DEFINE_THIS
iflags=$IFLAGS     # Location of include files - shouldn't be needed
f77=$F77           # Fortran 77 serial compiler (must be in execution path)
f90=$F90           # Fortran 90 serial compiler (must be in execution path)
mpif90=$MPIF90     # Fortran 90 parallel compiler (must be in execution path)
fflags=$FFLAGS     # Flags for Fortran 77 and 90 compilers
fflags_nomain=$FFLAGS_NOMAIN # Flags for linking Fortran sources with main in a different language
fflags_noopt=$FFLAGS_NOOPT # as FFLAGS With optimization disabled
f90flags=$F90FLAGS # Flags for Fortran 90 compiler only
ld=$LD             # Loader (must be in the execution path)
ldflags=$LDFLAGS   # Flags for loader
ld_libs=$LD_LIBS   # Additional libraries
blas_libs=$BLAS_LIBS     # blas library - specify e.g. /my/blas/lib/libmyblas.a
                         # or -L/my/blas/lib -lmyblas
lapack_libs=$LAPACK_LIBS # lapack library, similar to above
fft_libs=$FFT_LIBS       # FFT libraries - may depend upon DFLAGS
mpi_libs=$MPI_LIBS       # MPI libraries - shouldn't be needed
mass_libs=$MASS_LIBS     # MASS libraries (IBM only)
libdirs=$LIBDIRS         # Where to look for libraries (e.g. /my/blas/lib)
scalapack_libs=$SCALAPACK_LIBS # scalapack libs
scalapack_dir=$SCALAPACK_LIB  # Where to look for scalapack libs
blacs_dir=$BLACS_LIB          # Where to look for libblacs.a
ar=$AR                   # ar (shouldn't be needed)
arflags=$ARFLAGS         # Flags for ar (as above)
extlib_flags=$EXTLIB_FLAGS # Flags for internal copies of lapack and blas
# configure for current directory by default
if test "$topdir" = "" ; then topdir="`pwd`" ; fi

# check system type (no cross-compilation for now)
AC_CANONICAL_BUILD

# many HPC systems are configured so that running parallel programs
# interactively is disabled: on those systems, AC_PROG_FC / _F77 / _CC
# would fail because they can't run the compiled executables.
# to work around that, let's pretend we are cross-compiling even if we aren't
# !!! this relies on undocumented Autoconf behavior !!!

# This is used to distinguish between true and fake cross compilation
# (only on NEC SX8 actually)
if test "$host" != "" ; then ranlib=echo; fi

cross_compiling=yes
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
                i*86-apple-darwin* )    arch=mac686 ;;
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
echo detected architecture... $arch

# candidate archiver and archiver flags
try_ar="ar"
try_arflags="ruv"
try_iflags="-I../include"

# candidate fortran compilers good for all cases
try_mpif90="mpif90"
try_f90="gfortran g95 f90"

# add all needed -D options to try_dflags
try_dflags=""

# candidate compilers and flags based on architecture
case $arch in
ia32 | ia64 | x86_64 )
        try_f90="ifort pgf90 pathf95 sunf95 openf95 $try_f90"
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
        try_arflags="-X64 ruv"
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
        sxopt=`echo $host|awk '{print substr($1,1,3)}'`
        echo $sxopt $host
        try_mpif90="sxmpif90"
        try_f90="sxf90"
        try_ar="sxar"
        try_arflags="rv"
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
        try_arflags="ruv"
        ;;
# PowerPC MareNostrum
ppc64-mn )
        try_f90="xlf90_r"
        try_dflags="-D__XLF"
        try_arflags="ruv"
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

if test "$enable_environment" -eq 1 ;
then
        try_dflags="$try_dflags -D__ENVIRONMENT"
fi

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
        eko_version=`$mpif90 -v 2>&1 | grep "EKOPath"` # obsolescent
        pathf95_version=`$mpif90 -v 2>&1 | grep "PathScale(TM)"`
        gfortran_version=`$mpif90 -v 2>&1 | grep "gcc version"`
        #
        if test "$ifort_version" != ""
        then
                version=`$mpif90 -V 2>&1 | grep Version |
                         sed 's/.*Version//' | awk '{print $1}'`
                ifort_version=`echo $version | sed 's/\..*//'`
                echo "${ECHO_T}ifort $version"
                f90_in_mpif90="ifort"
                if test "$ifort_version" -gt 8; then
                # flags for MKL - ifort 9 and later
                   MKL_LIBS=""
                   if test "$ifort_version" -gt 9; then
                        MKL_FLAGS="-static-intel"
                   else
                        MKL_FLAGS="-i-static"
                   fi
                else
                # flags for MKL - ifort 8 and earlier, obsolescent
                   MKL_LIBS="-lguide -lpthread"
                   MKL_FLAGS=""
                fi
        elif test "$sunf95_version" != ""
        then
                version=`echo $sunf95_version | awk '{print $5}'`
                echo "${ECHO_T}sunf95 $version"
                f90_in_mpif90="sunf95"
        elif test "$openf95_version" != ""
        then
                version=`echo $openf95_version | awk '{print $5}'`
                echo "${ECHO_T}openf95 $version"
                f90_in_mpif90="openf95"
        elif test "$pgf_version" != ""
        then
                version=`echo $pgf_version | awk '{print $2}'`
                echo "${ECHO_T}pgf90 $version"
                f90_in_mpif90="pgf90"
        elif test "$eko_version" != ""
        then
# obsolescent
                version=`echo $eko_version | awk '{print $6}'`
                echo "${ECHO_T}pathf95 $version"
                f90_in_mpif90="pathf95"
        elif test "$g95_version" != ""
        then
                version=`echo $g95_version | awk '{print $3}'`
                echo "${ECHO_T}g95 $version"
                f90_in_mpif90="g95"
        elif test "$pathf95_version" != ""
        then
                version=`echo $pathf95_version | awk '{print $5}'`
                echo "${ECHO_T}pathf95 $version"
                f90_in_mpif90="pathf95"
        elif test "$gfortran_version" != ""
        then
                version=`echo $gfortran_version | awk '{print $3}'`
                echo "${ECHO_T}gfortran $version"
                f90_in_mpif90="gfortran"
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
        f90_version=ifort
    elif $f90 -V 2>&1 | grep -q "Sun Fortran" ; then
        f90_version=sunf95
    elif $f90 -V 2>&1 | grep -q "^Open64" ; then
        f90_version=openf95
    elif $f90 -V 2>&1 | grep -q "^pgf" ; then
        f90_version=pgf
    elif $f90 -v 2>&1 | grep -q "EKOPath" ; then # obsolescent
        f90_version=pathf95
    elif $f90 -v 2>&1 | grep -q "g95" ; then
        f90_version=g95
    elif $f90 -version 2>&1 | grep -q "PathScale(TM)" ; then
        f90_version=pathf95
    elif $f90 -v 2>&1 | grep -q "gcc version" ; then
        f90_version=gfortran
    elif $f90 -V 2>&1 | grep -q "Cray Fortran" ; then
        f90_version=crayftn
    else
        echo $ECHO_N "unknown, leaving as... $ECHO_C"
        f90_version=$f90
    fi
    echo $f90_version
    ;;
* )
    f90_version=$f90
    ;;
esac



# candidate C and f77 compilers good for all cases
try_cc="cc gcc"
try_f77="$f90"

case "$arch:$f90_version" in
*:ifort* )
        try_cc="icc ecc $try_cc"
        ;;
*:pgf90 )
        try_cc="pgcc $try_cc"
        try_f77="pgf77 $f90"
        ;;
*:pathf95 )
        try_cc="pathcc $try_cc"
        ;;
*:sunf95 )
        try_cc="suncc $try_cc"
        ;;
*:openf95 )
        try_cc="opencc $try_cc"
        ;;
aix:*xlf*_r )
        try_cc="xlc_r $try_cc"
        try_f77="xlf_r $f90"
        ;;
aix:*xlf* )
        try_cc="xlc $try_cc"
        try_f77="xlf $f90"
        ;;
sparc:* | solaris:sunf95 )
        try_f77="f77 $f90"
        ;;
cray*:* )
        try_cc="cc"
        ;;
necsx:* )
        try_cc="sxcc"
        ;;
ppc64-bg*:*xlf90_r )
        try_cc="bgxlc_r"
        try_f77="bgxlf_r"
        ;;
ppc64-bg*:*xlf90 )
        try_cc="bgxlc"
        try_f77="bgxlf"
        ;;
ppc64:*xlf* | ppc64-mn:*xlf* )
        try_cc="xlc_r $try_cc"
        try_f77="xlf_r $f90"
        ;;
esac

# check serial C compiler
if test "$env_cc" = "" ; then cc="$try_cc" ; else cc="$env_cc"; fi
AC_PROG_CC($cc)
cc=$CC

echo setting CC... $cc

# check size of pointers to int - needed to decide the size of integer
# arrays in fortran holding C pointers for FFTW

AC_CHECK_SIZEOF([int *])
SIZEOF_INT_P=$ac_cv_sizeof_int_p
AC_SUBST(SIZEOF_INT_P)
AC_CONFIG_FILES([include/fft_defs.h:include/fft_defs.h.in])

# check if the structure mallinfo is present in malloc.h
AC_CHECK_HEADER(malloc.h,have_malloc_h=1,have_malloc_h=0, )
if test "$have_malloc_h" -ne 0
then
AC_CHECK_MEMBER([struct mallinfo.arena],
                [AC_DEFINE(HAVE_MALLINFO)],
                ,
                [#include <malloc.h>])

fi

# check serial Fortran 77 compiler (use F77 if it was set)
if test "$f77" = "" ; then f77="$try_f77" ; fi
AC_PROG_F77($f77)
f77=$F77

echo setting F77... $f77

# check Fortran compiler flags
# have_cpp=0: use external C preprocessing for fortran code
# have_cpp=1: use C-like preprocessing in fortran compiler
have_cpp=1
xlf_flags=0

echo using F90... $f90

case "$arch:$f90_version" in
ia32:ifort* | ia64:ifort* | x86_64:ifort* | mac686:ifort* | crayxt*:ifort* )
        try_fflags="-O2 -assume byterecl -g -traceback -par-report0 -vec-report0"
        if test "$use_debug" -eq 1; then
            try_fflags="$try_fflags -fpe0 -CB"
        fi
	    try_fflags_nomain="-nofor_main"
        try_fflags_openmp="-openmp"
        try_f90flags="\$(FFLAGS) -nomodule"
        try_fflags_noopt="-O0 -assume byterecl -g -traceback"
        try_ldflags=""
        try_ldflags_static="-static"
        try_ldflags_openmp="-openmp"
        try_dflags="$try_dflags -D__INTEL"
        pre_fdflags="-fpp "
        ;;
ia32:pgf* | ia64:pgf* | x86_64:pgf* )
	    try_fflags_nomain="-Mnomain"
        try_fflags="-fast -r8"
        try_fflags_openmp="-mp"
        try_f90flags="-fast -r8 -Mcache_align"
        try_fflags_noopt="-O0"
        try_ldflags=""
        try_ldflags_openmp="-mp"
        try_ldflags_static="-Bstatic"
        try_dflags="$try_dflags -D__PGI"
        have_cpp=0
        ;;
ia32:path* | ia64:path* | x86_64:path* )
        try_fflags="-march=auto -O2"
        try_f90flags="\$(FFLAGS)"
        try_fflags_noopt="-O0"
        try_ldflags=""
        try_ldflags_static="-static"
        have_cpp=0
        ;;
*:g95 )
        if test "$use_debug" -eq 1; then
            try_fflags="-O3 -g -freal=nan -finteger=12345678 -flogical=none -cpp"
        else
            try_fflags="-O3 -cpp"
        fi
        try_f90flags="\$(FFLAGS)"
        try_fflags_noopt="-O0 -cpp"
        try_ldflags=""
        try_ldflags_static="-static"
        ;;
*:*gfortran )
        if test "$use_debug" -eq 1; then
            try_fflags="-O3 -g  -Wall -fbounds-check -frange-check"
        else
            try_fflags="-O3 -g"
        fi
        try_fflags_openmp="-fopenmp"
        try_f90flags="\$(FFLAGS) -x f95-cpp-input"
        try_fflags_noopt="-O0 -g"
        try_ldflags="-g -pthread"
        try_ldflags_openmp="-fopenmp"
        try_dflags="$try_dflags -D__GFORTRAN -D__STD_F95"
        try_ldflags_static="-static"
        ;;
*:sunf95 )
        try_fflags="-O4"
        try_fflags_openmp="-openmp"
        try_f90flags="\$(FFLAGS) -fpp"
        try_fflags_noopt="-O0"
        try_ldflags="-fast"
        try_ldflags_static="-Bstatic"
        imod="-M"
        ;;
*:openf95 )
        try_fflags="-O3"
        try_f90flags="\$(FFLAGS) -ftpp"
        try_fflags_noopt="-O0"
        try_ldflags=""
        imod="-I"
        ;;
aix:*xlf* )
        if test "$use_debug" -eq 1; then
            try_fflags="-q64 -qalias=noaryovrlp -g -C \
-qarch=auto -qtune=auto -qdpc -Q -qalias=nointptr"
        else
            try_fflags="-q64 -qalias=noaryovrlp -O3 -qstrict \
-qarch=auto -qtune=auto -qdpc -Q -qalias=nointptr"
        fi
        try_fflags_openmp="-qsmp=omp"
        try_f90flags="\$(FFLAGS) -qsuffix=cpp=f90 -qfree=f90"
        try_fflags_noopt="-q64 -O0"
        try_ldflags="-q64"
        try_ldflags_openmp="-qsmp=omp"
        # try_ldflags_static="-bstatic"
        pre_fdflags="-WF,"
        xlf_flags=1
        ;;
solaris:sunf95 )
        try_fflags="-fast -O2 -fpp"
        try_f90flags="\$(FFLAGS)"
        try_fflags_noopt="-O0 "
        try_ldflags=""
        imod="-M"
        ;;
sparc:f90 )
        try_fflags="-fast -O1 -nodepend -xvector=no -xchip=ultra3 \
-xarch=v8plusb -xlic_lib=sunperf"
        try_f90flags="\$(FFLAGS)"
        try_fflags_noopt="-O0 -xlic_lib=sunperf"
        try_ldflags=""
        imod="-M"
        have_cpp=0
        ;;
crayxt*:cray* )
        try_fflags_nomain=""
        #NOTE: by default OpenMP is always ON (see crayftn man page)
        try_fflags_openmp="-homp"
        try_fflags="-O2"
        #NOTE: add '-rm' to get messages from crayftn about why
        #      optimizations have not been applied
        try_f90flags="-O3,fp3 -f free"
        try_fflags_noopt="-O0"
        try_ldflags_openmp="-homp"
        try_ldflags="-v"
        try_ldflags_static="-static"
        try_dflags="$try_dflags -D__CRAY"
        have_cpp=0
        ;;
crayxt*:pgf* )
# see comment above for pgf*
	    try_fflags_nomain="-Mnomain"
        try_fflags_openmp="-mp"
        try_fflags="-O3 -r8"
        try_f90flags="-fast -Mcache_align -r8 -Mpreprocess"
        try_fflags_noopt="-O0"
        try_ldflags_openmp="-mp"
        try_ldflags="-v"
        try_dflags="$try_dflags -D__PGI -D__IOTK_WORKAROUND1"
        have_cpp=1
        ;;
crayxt*:pathf* )
        try_fflags="-march=auto -O2 -cpp"
        try_f90flags="\$(FFLAGS)"
        try_fflags_noopt="-O0"
        try_ldflags=""
        try_ldflags_static="-static"
        have_cpp=1
        ;;
necsx:* )
        try_fflags='      -float0 -Cvopt -eab -R5 -Wf,-Ncont,-A dbl4,-P nh,-ptr byte,-pvctl noifopt loopcnt=9999999 expand=12 fullmsg vwork=stack,-fusion,-O noif,-init stack=nan heap=nan'
        try_f90flags='  -f2003  -float0 -Cvopt -eab -R5 -Wf,-Ncont,-A dbl4,-P nh,-ptr byte,-pvctl noifopt loopcnt=9999999 expand=12 fullmsg vwork=stack,-fusion,-O noif,-init stack=nan heap=nan'
        try_f90flags="-$sxopt $try_f90flags"
        try_fflags_noopt='-float0   '
        try_f90flags_noopt='-f2003 -float0 -eab -R5 -C debug  -Wf,-Ncont,-A dbl4,-P nh ,ptr byte,-init stack=nan heap=nan'
        try_f90flags_noopt="$try_f90flags_noopt"
        try_f90flags_inline='-f2003  -float0 -Cvopt -eab -R5 -pi noauto incdir exp=w0gauss -Wf,-Ncont,-A dbl4,-P nh,-ptr byte,-pvctl noifopt loopcnt=9999999 expand=12 fullmsg vwork=stack,-fusion,-O noif,-init stack=nan heap=nan'
        try_f90flags_inline="$try_f90flags_inline"
        try_ldflags_static='-P static'
        try_ldflags='-Wl,-f zero'
        try_ldflags="-p $try_ldflags"
        pre_fdflags=""
        ;;

ppc64:*xlf* )
    if test "$use_debug" -eq 1; then
        try_fflags="-g -C -qsuffix=cpp=f90 -qdpc -qalias=nointptr -Q"
    else
        try_fflags="-q64 -qthreaded -O4 -qsuffix=cpp=f90 -qdpc -qalias=nointptr -Q"
    fi
        try_f90flags="\$(FFLAGS) -qfree=f90"
        try_fflags_noopt="-q64 -qthreaded -O0"
        try_ldflags="-q64 -qthreaded"
        pre_fdflags="-WF,"
        xlf_flags=1
        ;;
ppc64-mn:*xlf* )
    if test "$use_debug" -eq 1; then
        try_fflags="-g -C -q64 -qstrict -qsuffix=cpp=f90 -qdpc -qalias=nointptr -Q -qtune=ppc970 -qarch=ppc970 -qcache=auto -qhot=vector,simd -qenablevmx"
    else
        try_fflags="-O3 -q64 -qstrict -qsuffix=cpp=f90 -qdpc -qalias=nointptr -Q -qtune=ppc970 -qarch=ppc970 -qcache=auto -qhot=vector,simd -qenablevmx"
    fi
        try_f90flags="\$(FFLAGS) -qfree=f90"
        try_fflags_noopt="-O0 -q64"
        try_ldflags=""
        pre_fdflags="-WF,"
        xlf_flags=1
        ;;
ppc64-bg:*xlf* )
    if test "$use_debug" -eq 1; then
        try_fflags="-q32 -qalias=noaryovrlp:nointptr -g -C -qdpc=e"
    else
        try_fflags="-q32 -qalias=noaryovrlp:nointptr -O3 -qstrict -qdpc=e"
    fi
        try_fflags_openmp="-qsmp=omp -qthreaded"
        try_f90flags="\$(FFLAGS) -qsuffix=cpp=f90"
        try_fflags_noopt="-q32 -O0"
        try_ldflags="-q32"
        try_ldflags_openmp="-qsmp=omp -qthreaded"
        pre_fdflags="-WF,"
        xlf_flags=1
        ;;
ppc64-bgq:*xlf* )
    if test "$use_debug" -eq 1; then
        try_fflags="-qalias=noaryovrlp:nointptr -g -C -qdpc=e"
    else
        try_fflags="-qalias=noaryovrlp:nointptr -O3 -qstrict -qdpc=e -qarch=qp -qtune=qp"
    fi
        try_fflags_openmp="-qsmp=noauto:omp -qtm -qthreaded"
        try_f90flags="\$(FFLAGS) -qsuffix=cpp=f90"
        try_fflags_noopt="-O0"
        try_ldflags=""
        try_ldflags_openmp="-qstatic -qsmp=noauto:omp -qtm -qthreaded"
        pre_fdflags="-WF,"
        xlf_flags=1
        ;;

* )
        # unknown, try these
        try_fflags="-O1"
        try_f90flags="\$(FFLAGS)"
        try_fflags_noopt="-O0"
        try_ldflags=""
        have_cpp=0
        ;;

esac

if test "$use_shared" -eq 0 ; then
  try_ldflags="$try_ldflags $try_ldflags_static" ; fi

if test "$use_openmp" -eq 1 ; then
  try_f90flags="$try_f90flags $try_fflags_openmp"
  try_fflags="$try_fflags $try_fflags_openmp"
  try_ldflags="$try_ldflags $try_ldflags_openmp"
fi

if test "$fflags" = ""   ; then fflags=$try_fflags     ; fi
if test "$f90flags" = "" ; then f90flags=$try_f90flags ; fi
if test "$fflags_noopt" = ""   ; then fflags_noopt=$try_fflags_noopt     ; fi
if test "$fflags_nomain" = ""   ; then fflags_nomain=$try_fflags_nomain     ; fi

echo setting FFLAGS... $fflags
echo setting F90FLAGS... $f90flags
echo setting FFLAGS_NOOPT... $fflags_noopt
if test "$fflags_nomain" != "" ; then echo setting FFLAGS_NOMAIN... $fflags_nomain ; fi

if test "$imod" = "" ; then imod="-I" ; fi

# tentative C and loader flags, good for many cases
try_cflags="-O3"
c_ldflags=""
try_cpp="cpp"

case "$arch:$cc" in
*:pgcc )
        try_cflags="-fast"
        ;;
# pathcc -E seems to give problems when preprocessing iotk
#*:pathcc )
#       try_cpp="pathcc -E"
#       ;;
aix:xlc* | aix:cc )
        try_cflags="-q64 -O2"
        c_ldflags="-q64"
        ;;
*:suncc  )
        try_cflags="-fast -O"
        ;;
sparc:cc )
        try_cflags="-fast -dalign -xchip=ultra3 -xarch=v8plusb \
-xlic_lib=sunperf"
        try_cpp="fpp"
        ;;
crayxt*:cc )
        # Actually we need something like is done for ftn to detect 
        # the proper compiler used (NdFilippo)
        try_cflags="-O3"
        ;;
necsx:* )
        #try_cflags="-D__SX6 \$(IFLAGS) \$(MODFLAGS)"
        try_cflags=""
        try_ar="sxar"
        ;;
ppc64-mn:* )
        try_cflags="-O3 -q64"
        ;;
ppc64-bg:* )
        try_cflags="-O3 -q32"
        ;;
ppc64-bgq:* )
        try_cflags="-O3"
        ;;
ppc64:xlc*)
        try_cflags="-O3 -q64 -qthreaded"
        c_ldflags="-q64"
        ;;

esac
if test "$cflags" = "" ; then cflags=$try_cflags ; fi
echo setting CFLAGS... $cflags

# preprocessor - try cpp in all cases; the preprocessor returned by
# AC_PROG_CPP
# may sometimes refuse to preprocess fortran files
if test "$cpp" = "" ; then cpp=$try_cpp; fi
# if test "$cpp" = "" ; then cpp=$CPP;     fi
echo setting CPP... $cpp

echo $ECHO_N "setting CPPFLAGS... $ECHO_C"
case $cpp in
        cpp)  try_cppflags="-P -C -traditional" ;;
        fpp)  try_cppflags="-P -C "              ;;
        *)    try_cppflags=""                ;;
esac
if test "$cppflags" = "" ; then cppflags=$try_cppflags ; fi
echo "${ECHO_T}$cppflags"

# linker and archiver
# note that from this point on, further additions to
# linker flags should be added to ldflags rather than try_ldflags
if test "$ld" = ""       ; then ld="$mpif90"           ; fi
if test "$ldflags" = ""  ; then ldflags="$try_ldflags" ; fi
if test "$ar" = ""       ; then ar="$try_ar"           ; fi
if test "$arflags" = ""  ; then arflags="$try_arflags" ; fi
echo setting LD... $ld
echo setting LDFLAGS... $ldflags
echo setting AR... $ar
echo setting ARFLAGS... $arflags

# compilation rules
AC_PROG_MAKE_SET
echo $ECHO_N "checking whether Fortran files must be preprocessed... $ECHO_C"
if test "$have_cpp" -ne 0
then
        f90rule="\$(MPIF90) \$(F90FLAGS) -c \$<"
        echo "${ECHO_T}no"
else
        f90rule="\$(CPP) \$(CPPFLAGS) \$< -o \$*.F90 ; \\
        \$(MPIF90) \$(F90FLAGS) -c \$*.F90 -o \$*.o"
        echo "${ECHO_T}yes"
fi

# compilation flags for all subsequent tests
# remove all $(...) because at least one compiler doesn't like them
# but if f90flags contains $(FFLAGS), substitute it
test_cflags="`echo $cflags | sed 's/\$([[^)]]*)//g'`"
test_cppflags="$test_cflags"
if test "`echo $f90flags | grep '$(FFLAGS)'`" != ""
then
        test_fflags="`echo $fflags $f90flags | sed 's/\$([[^)]]*)//g'`"
else
        test_fflags="`echo $f90flags | sed 's/\$([[^)]]*)//g'`"
fi
test_ldflags="`echo $ldflags | sed 's/\$([[^)]]*)//g'`"

# find Fortran to C wrappers

if test "$check_wrappers" -ne 0; then
   AC_F77_WRAPPERS
fi

AC_LANG_PUSH(Fortran 77)
F77=$f90 # use Fortran 90 actually

# Search for libraries
# Flags for needed libraries
have_blas=0
have_lapack=0
have_fft=0
have_mpi=0
have_scalapack=0
# Flags for machine-specific libraries
have_acml=0
have_atlas=0
have_essl=0
have_mkl=0


# check for blas
# supported vendor replacements:
#   mkl and acml on Intel/AMD architectures
#   essl on aix
#   SUNperf on sparc
# atlas is used over blas if available
# internal version is used if none is found

if test "$use_internal_blas" -eq 0
then
   if test "$blas_libs" = ""
   then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        case "$arch:$f90" in

        x86_64:pgf* | x86_64:path* | x86_64:openf95 | crayxt*:* )
                # check for acml - note that it contains lapack as well
                try_libdirs="/opt/acml*/pathscale64/lib/"
                try_libdirs="$ld_library_path $libdirs $try_libdirs"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_dgemm # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        
                        if test "$use_openmp" -eq 0; then
                                AC_SEARCH_LIBS(dgemm, acml, have_blas=1 have_lapack=1
                                    have_acml=1 blas_libs="$try_loption $LIBS")
                        else
                                AC_SEARCH_LIBS(dgemm, acml_mp, have_blas=1 have_lapack=1
                                    have_acml=1 blas_libs="$try_loption $LIBS")                        
                        fi
                        
                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi
                done
                ;;

        ia64:* )
                # check for mkl (in several directories)
                try_libdirs="/opt/intel/Compiler/*/*/mkl/lib/64
                             /opt/intel/mkl/*/lib/64
                             /opt/intel/mkl*/lib/64
                             /opt/intel/mkl/lib"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_dgemm # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=" "
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$MKL_FLAGS $test_ldflags $try_loption"
                        LIBS="$MKL_LIBS"
                        #
                        # should work for recent MKL versions only
                        #
                        if test "$use_openmp" -eq 0; then
                           if test "$f90" = "g95" -o "$f90" = "gfortran" ; then
     			      AC_SEARCH_LIBS(dgemm, mkl_gf_ipf, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_sequential -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   else
     			      AC_SEARCH_LIBS(dgemm, mkl_intel_ipf, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_sequential -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   fi
                        else
                           if test "$f90" = "g95" -o"$f90" = "gfortran"; then
     			      AC_SEARCH_LIBS(dgemm, mkl_gf_ipf, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_gnu_thread -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   else
     			      AC_SEARCH_LIBS(dgemm, mkl_intel_ipf, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_intel_thread -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   fi
                        fi
                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi
                done
                ;;

        *:sunf95 )
                # check for acml - note that it contains lapack as well
                if test "$arch" = "x86_64"
                then
                        try_libdirs="/usr/local/sunstudio*/lib/amd64/"
                else
                        try_libdirs="/usr/local/sunstudio*/lib/"
                fi
                try_libdirs="$libdirs $ld_library_path $try_libdirs"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_dgemm # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        AC_SEARCH_LIBS(dgemm, sunperf, have_blas=1 have_lapack=1
                                blas_libs="$try_loption $LIBS")
                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi
                done
                ;;

        x86_64:* )
                try_libdirs="/opt/intel/composer*/mkl/lib/intel64
                             /opt/intel/Compiler/*/*/mkl/lib/em64t
                             /opt/intel/mkl/*/lib/em64t
                             /opt/intel/mkl*/lib/em64t
                             /opt/intel/mkl/lib"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_dgemm # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=" "
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$MKL_FLAGS $test_ldflags $try_loption"
                        LIBS="$MKL_LIBS"
                        #
                        # should work for recent MKL versions only
                        #
                        if test "$use_openmp" -eq 0; then
                           if test "$f90" = "g95" -o "$f90" = "gfortran" ; then
     			      AC_SEARCH_LIBS(dgemm, mkl_gf_lp64, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_sequential -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   else
     			      AC_SEARCH_LIBS(dgemm, mkl_intel_lp64, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_sequential -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   fi
                        else
                           if test "$f90" = "g95" -o "$f90" = "gfortran" ; then
     			      AC_SEARCH_LIBS(dgemm, mkl_gf_lp64, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_gnu_thread -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   else
     			      AC_SEARCH_LIBS(dgemm, mkl_intel_lp64, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_intel_thread -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   fi
                        fi
                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi
                done
                ;;

        ia32:* )
                # check for mkl (in several directories)
                try_libdirs="/opt/intel/composer*/mkl/lib/ia32
                             /opt/intel/Compiler/*/*/mkl/lib/32
                             /opt/intel/mkl/*/lib/32
                             /opt/intel/mkl*/lib/32
                             /opt/intel/mkl/lib"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_dgemm # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption="-L "
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$MKL_FLAGS $test_ldflags $try_loption"
                        LIBS="$MKL_LIBS"
                        #
                        # should work for recent MKL versions only
                        #
                        if test "$use_openmp" -eq 0; then
                           if test "$f90" = "g95" -o "$f90" = "gfortran"; then
     			      AC_SEARCH_LIBS(dgemm, mkl_gf, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_sequential -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   else
     			      AC_SEARCH_LIBS(dgemm, mkl_intel, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_sequential -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   fi
                        else
                           if test "$f90" = "g95" -o "$f90" = "gfortran" ; then
     			      AC_SEARCH_LIBS(dgemm, mkl_gf, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_gnu_thread -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   else
     			      AC_SEARCH_LIBS(dgemm, mkl_intel, 
                                 have_blas=1 have_mkl=1 
                                 blas_libs="$try_loption $LIBS -lmkl_intel_thread -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core)
			   fi
                        fi
                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi

                done
                ;;

        aix:* )
                # check for essl
                unset ac_cv_search_dgemm # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS=""
                AC_SEARCH_LIBS(dgemm, essl, have_blas=1
                               blas_libs="$LIBS" )
                # notice that some IBM machines may not need -lessl
                # to load blas so the above test may fail
                if test "`echo $blas_libs | grep essl`" != ""
                then
                    have_essl=1
                    try_dflags="$try_dflags -D__ESSL"
                fi
		# we need esslsmp for hybrid (MPI+OpenMP) build
		if test "$have_essl"="1"; then
		    if test "$use_openmp" -ne 0 ; then
		         blas_libs="-lesslsmp"
		    fi
		fi
                ;;

        sparc:* | solaris:* )
                # check for SUNperf library
                unset ac_cv_search_dgemm # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS=""
                AC_SEARCH_LIBS(dgemm, sunperf, have_blas=1 have_lapack=1
                               blas_libs="-xlic_lib=sunperf $LIBS")
                ;;
        necsx:* )
                #sx5-nec or sx6-nec or sx8-nec: check in (/SX)/usr/lib
                #sx8-nec-idris: check in /SX/opt/mathkeisan/inst/lib0
                try_libdirs="/SX/usr/lib /SX/opt/mathkeisan/inst/lib0"
                for dir in none $try_libdirs
                do
                        unset ac_cv_search_dgemm # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        AC_SEARCH_LIBS(dgemm, blas, have_blas=1
                                       blas_libs="$try_loption $LIBS")
                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi
                 done
                 ;;
        ppc64:* )
                # check for essl
                unset ac_cv_search_dgemm # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS=""
                AC_SEARCH_LIBS(dgemm, essl, have_blas=1
                               blas_libs="$LIBS" )
                # notice that some IBM machines may not need -lessl
                # to load blas so the above test may fail
                if test "`echo $blas_libs | grep essl`" != ""
                then
                    have_essl=1
                    try_dflags="$try_dflags -D__LINUX_ESSL"
                fi
                # OBM:Yet another work-around if the above search 
                # returns "none required" 
                if test "$ac_cv_search_dgemm" = "none required"
                then
                    echo "There is no need for -lessl in this machine"
                    have_essl=1
                    try_dflags="$try_dflags -D__LINUX_ESSL"
                fi
		# we need esslsmp for hybrid (MPI+OpenMP) build
		if test "$have_essl"="1"; then
		    if test "$use_openmp" -ne 0 ; then
		         blas_libs="-lesslsmp"
		    fi
		fi
                ;;

        ppc64-*:* )
                # assume essl
                unset ac_cv_search_dgemm # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                have_blas=1
                have_essl=1
		# BlueGene: for some obscure reason there is no need to
		# specify a library path to have essl linked, while
		# in reality it is needed to specify where essl are
		if test "$arch"="ppc64-bg"; then
                   try_dflags="$try_dflags -D__LINUX_ESSL"
		   if test "$blas_libs"=""; then
		      if test "$use_openmp" -eq 0 ; then
		         blas_libs="-L/opt/ibmmath/essl/4.4/lib/ -lesslbg"
		      else
		         blas_libs="-L/opt/ibmmath/essl/4.4/lib/ -lesslsmpbg"
		      fi
		   fi
                else
                   try_dflags="$try_dflags -D__LINUX_ESSL"
		fi
                ;;
	mac686:ifort* )
                #This solution is tested with MacOs 10.6 and Intel 11.1
                #..and now MacOs 10.8.3 and Intel 13 
                try_libdirs="/Developer/opt/intel/Compiler/*/*/Frameworks/mkl/lib/universal
                             /opt/intel/Compiler/*/*/Frameworks/mkl/lib/universal
                             /opt/intel/mkl*/lib/em64t
                             /opt/intel/mkl/lib"
		try_libdirs="$libdirs $try_libdirs $ld_library_path"

		for dir in none $try_libdirs
		do
			unset ac_cv_search_dgemm # clear cached value
			if test "$dir" = "none"
			then
				try_loption=""
			else
				echo $ECHO_N "in $dir: " $ECHO_C
				try_loption="-L$dir"
			fi
			FFLAGS="$test_fflags"
			LDFLAGS="$MKL_FLAGS $test_ldflags $try_loption"
			LIBS="$MKL_LIBS"
                        # First, a by-the-apple-book search of MKL... >10.2 requires multiple libraries
                        # 64 bit is buggy as of 11.1.088
                        if test "$use_openmp" -eq 0; then
                        AC_SEARCH_LIBS(dgemm, mkl_intel,
                                       have_blas=1 have_mkl=1
                                       blas_libs="$try_loption $LIBS -lmkl_sequential -lmkl_core -lpthread"
                                       ldflags="$MKL_FLAGS $ldflags",
                                       echo "MKL not found",
                                       -lmkl_sequential -lmkl_core -lpthread)
			else 
                        AC_SEARCH_LIBS(dgemm, mkl_intel,
                                       have_blas=1 have_mkl=1
                                       blas_libs="$try_loption $LIBS -lmkl_intel_thread -lmkl_core -openmp -lpthread"
                                       ldflags="$MKL_FLAGS $ldflags",
                                       echo "MKL not found",
                                       -lmkl_intel_thread -lmkl_core -openmp -lpthread)
			fi
                        # 32 bit
                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi
      		done
		;;


        esac
        # blas not (yet) found: look for more possibilities
        if test "$have_blas" -eq 0
        then
        case "$f90" in
                pgf* )
                # check for PGI blas
                unset ac_cv_search_dgemm # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS=""
                AC_SEARCH_LIBS(dgemm, blas, have_blas=1 blas_libs="$LIBS")
                ;;
        esac
        fi

        if test "$have_blas" -eq 0
        then
                # check for atlas (in several directories)
                try_libdirs="/usr/local/lib"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_dgemm # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS="-latlas"
                        AC_SEARCH_LIBS(dgemm, f77blas, have_blas=1 have_atlas=1
                                       blas_libs="$try_loption $LIBS", , -lg2c)
                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi
                done
        fi

        # blas still not found

        if test "$have_blas" -eq 0
        then
                # check for blas (in several directories)
                try_libdirs="/usr/local/lib"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_dgemm # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        AC_SEARCH_LIBS(dgemm, blas-3 openblas blas, have_blas=1
                                       blas_libs="$try_loption $LIBS")
                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi
                done
        fi
   else
        # blas provided in BLAS_LIBS - not checked!
        have_blas=1
   fi
fi

# check for lapack
#
# same supported vendor replacements as for blas
# internal version is used if none is found
if test "$use_internal_lapack" -eq 0
   then
   if test "$lapack_libs" = ""
   then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        case "$arch:$f90" in

        ia32:* | ia64:*| x86_64:* )
                # check for mkl_lapack (if mkl found and acml not found)
                # OBSOLESCENT - recent versions of mkl contain lapack
                if test "$have_mkl" -ne 0 && test "$have_acml" -eq 0
                then
                        unset ac_cv_search_dspev ac_lib # clear cached value
                        FFLAGS="$test_fflags"
                        LIBS=" $blas_libs"
                        LDFLAGS="$MKL_FLAGS $test_ldflags"
                        AC_SEARCH_LIBS(dspev, mkl_lapack, have_lapack=1)
                        if test "$ac_lib" != "" ; then lapack_libs="-l$ac_lib"; fi
                fi
                ;;

        sparc:* )
                # check for SUNperf library
                unset ac_cv_search_dspev # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$blas_libs"
                AC_SEARCH_LIBS(dspev, sunperf, have_lapack=1
                               lapack_libs="-xlic_lib=sunperf $LIBS")
                ;;
        aix:* )
                # check for essl
                unset ac_cv_search_dspev # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$blas_libs"
                AC_SEARCH_LIBS(dspev, essl, have_lapack=1
                                lapack_libs="$try_loption $LIBS"
                                try_dflags="$try_dflags -D__ESSL")
                # essl may not have been found in previous test on blas
                if test "$have_lapack" -eq 1; then have_essl=1; fi
                ;;
        ppc64:* )
                # check for essl
                unset ac_cv_search_dspev # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$blas_libs"
                AC_SEARCH_LIBS(dspev, essl, have_lapack=1
                                lapack_libs="$try_loption $LIBS"
                                try_dflags="$try_dflags -D__LINUX_ESSL")
                # essl may not have been found in previous test on blas
                if test "$have_lapack" -eq 1; then have_essl=1; fi
                ;;

        necsx:* )
                #sx5-nec or sx6-nec or sx8-nec: check in (/SX)/usr/lib
                #sx8-nec-idris: check in /SX/opt/mathkeisan/inst/lib0
                try_libdirs="/SX/usr/lib /SX/opt/mathkeisan/inst/lib0"
                for dir in none $try_libdirs
                do
                        unset ac_cv_search_dspev # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption $blas_libs"
                        LIBS=""
                        AC_SEARCH_LIBS(dspev, lapack, have_lapack=1
                                       lapack_libs="$try_loption $LIBS")
                        if test "$ac_cv_search_dspev" != "no"
                        then break ; fi
                done
                ;;
        esac

        if test "$have_lapack" -eq 0
        then
                # check for lapack (in several directories)
                try_libdirs="/usr/local/lib"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_dspev # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS="$blas_libs"
                        AC_SEARCH_LIBS(dspev, lapack-3 lapack, have_lapack=1
                                       lapack_libs="$try_loption $LIBS")
                        if test "$ac_cv_search_dspev" != "no"; then
                           # essl must precede lapack (if present)
                           if test "$have_essl" -ne 0 ; then
                                lapack_libs="$blas_libs $lapack_libs"
                           fi
                           break
                        fi
                done
        fi

   else
        # lapack provided in LAPACK_LIBS - not checked!
        have_lapack=1
   fi
fi

# no blas library found, or internal blas required: use the built-in blas
# (blas_libs is used in the above lapack tests: do not move the following
# settings above lapack tests, which would seem a more logical place)

if test "$have_blas" -eq 0 -o "$use_internal_blas" -eq 1 ; then
    blas_libs="$topdir/BLAS/blas.a"
    blas_libs_switch="internal"
else
    blas_libs_switch="external"
fi

# no lapack library found, or incomplete lapack found (atlas, essl),
# or internal lapack esplicitly required

if test "$have_lapack" -eq 0 -o "$use_internal_lapack" -eq 1 ; then
    lapack_libs="$topdir/lapack-3.2/lapack.a"
    lapack_libs_switch="internal"
else
    if test "$have_essl" -eq 1 -o "$have_atlas" -eq 1 ; then
    # IBM essl or atlas: add missing lapack routines - must be loaded after lib
    # atlas: add missing lapack routines so as to complete atlas
    # note that some compilers do not like to have multiple symbols
      lapack_libs="$lapack_libs $topdir/lapack-3.2/lapack.a"
      lapack_libs_switch="internal"
    else
      lapack_libs_switch="external"
    fi
fi

echo setting BLAS_LIBS... $blas_libs
echo setting LAPACK_LIBS... $lapack_libs

# check for FFT libraries (no check for explicit openmp)
# supported vendor replacements:
#   essl on aix and some IBM linux machines
#   SUNperf on sparc
#   ASL/Mathkeisan on Nec
#   acml on amd
if test "$fft_libs" = "" && test "$use_openmp" -eq 0
then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        case "$arch" in
        aix )
                # check for essl
                unset ac_cv_search_dcft # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$fft_libs"
                AC_SEARCH_LIBS(dcft, essl, have_fft=1 fft_libs="$LIBS")
            ;;
        ppc64 | ppc64-mn )
                # check for essl
                unset ac_cv_search_dcft # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$fft_libs"
                AC_SEARCH_LIBS(dcft, essl, have_fft=1 fft_libs="$LIBS")
            ;;
        ppc64-bg | ppc64-bgq )
                # check for esslbg
                unset ac_cv_search_dcft # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$fft_libs $blas_libs"
                AC_SEARCH_LIBS(dcft, esslbg, have_fft=1 fft_libs="$LIBS")
            ;;
        sparc )
                # check for SUNperf FFT library on Sun Sparcs
                # but not on solaris PC! it is slower than FFTW
                unset ac_cv_search_zfft3i # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$libs"
                AC_SEARCH_LIBS(zfft3i, sunperf, have_fft=1
                               try_dflags="$try_dflags -D__SUNPERF"
                               fft_libs="-xlic_lib=sunperf $LIBS")
                ;;
        necsx )
                if test "$use_fft_mathkeisan" -ne 0
                then
                   #sx5-nec or sx6-nec or sx8-nec: check in (/SX)/usr/lib
                   #sx8-nec-idris: check in /SX/opt/mathkeisan/inst/lib0
                   try_libdirs="/SX/usr/lib /SX/opt/mathkeisan/inst/lib0"
                   #check for Mathkeisan (Cray simple precision )
                   #search for initialization subroutine
                   echo $ECHO_N  "Searching in Mathkeisan" $ECHO_C
                   for dir in none $try_libdirs
                   do
                        unset ac_cv_search_zftfax # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        AC_SEARCH_LIBS(zftfax, fft, have_fft=1
                             try_dflags="$try_dflags try_dflags_fft_mathkeisan"
                                       fft_libs="$try_loption $LIBS")
                        if test "$ac_cv_search_zftfax" != "no"
                        then break ; fi
                   done
                fi
                if test "$use_fft_asl" -ne 0
                then
                   #check for asl in (/SX)/usr/lib
                   try_libdirs="/SX/usr/lib"
                   #search for initialization subroutine
                   echo $ECHO_N  "Searching in Asl" $ECHO_C
                   for dir in none $try_libdirs
                   do
                        unset ac_cv_search_zfc3cl # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        AC_SEARCH_LIBS(zfc3cl, asl, have_fft=1
                             asl_libs="$try_loption $LIBS"
                             try_dflags="$try_dflags $try_dflags_fft_asl"
                             fft_libs="$fft_libs $asl_libs")
                        if test "$ac_cv_search_zfc3cl" != "no"
                        then break ; fi
                   done
                fi
                if test "$use_fft_para" -ne 0
                then
                   try_dflags="$try_dflags $try_dflags_fft_para"
                fi
                ;;
        esac

fi

if test "$have_fft" -eq 0 && test "$use_openmp" -eq 0
then

        # check for fftw v.3 (in several directories)
        try_libdirs="/usr/local/lib"
        try_libdirs="$libdirs $try_libdirs $ld_library_path"

        for dir in none $try_libdirs
        do
                unset ac_cv_search_dfftw_execute_dft # clear cached value
                if test "$dir" = "none"
                then
                        try_loption=
                else
                        echo $ECHO_N "in $dir: " $ECHO_C
                        try_loption="-L$dir"
                fi
                CFLAGS="$test_cflags"
                CPPFLAGS="$test_cppflags"
                LDFLAGS="$c_ldflags $try_loption"
                LIBS="$fft_libs"
                AC_SEARCH_LIBS(dfftw_execute_dft, fftw3, have_fft=1
                               fft_libs="$try_loption $LIBS", , -lm)
                if test "$ac_cv_search_dfftw_execute_dft" != "no"
                then
                        try_dflags="$try_dflags -D__FFTW3"
                        break
                fi
        done

fi
echo setting FFT_LIBS... $fft_libs

# if no valid FFT library was found, use the local copy
# (This happens also if OpenMP is enabled...)
if test "$have_fft" -eq 0
then
        case "$arch" in
        ppc64-bg | ppc64-bgq )
            try_dflags="$try_dflags -D__LINUX_ESSL"
            ;;
        * )
		    try_dflags="$try_dflags -D__FFTW"
            ;;
        esac
fi

F77=$f90
FFLAGS="$test_fflags"
LDFLAGS="$test_ldflags"
# check for mass on aix
if test "$mass_libs" = ""
then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        case "$arch" in
        aix | ppc64-bg )
                # check for mass (in several directories)
                try_libdirs="/opt/ibmcmp/xlmass/bg/7.3/bglib64 /opt/ibmcmp/xlmass/bg/4.4/bglib /cineca/lib /cineca/lib/mass"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_vexp # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        AC_SEARCH_LIBS(vexp, massvp4 massv, , , -lmass)
                        if test "$ac_cv_search_vexp" = "-lmassvp4" \
                                -o "$ac_cv_search_vexp" = "-lmassv"
                        then mass_libs="$try_loption $ac_cv_search_vexp -lmass"
                        fi
                        if test "$ac_cv_search_vexp" != "no" ; then break ; fi
                done
                ;;
        ppc64-bgq )
                # check for mass (in several directories)
                try_libdirs="/opt/ibmcmp/xlmass/bg/7.3/bglib64"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_vexp # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        AC_SEARCH_LIBS(vexp, massv, , , -lmass_simd)
                        if test "$ac_cv_search_vexp" = "-lmassv"
                        then mass_libs="$try_loption $ac_cv_search_vexp -lmass_simd"
                        fi
                        if test "$ac_cv_search_vexp" != "no" ; then break ; fi
                done
                ;;

        ppc64* )
                # check for mass (in several directories)
                try_libdirs="/usr/local/lib /opt/ibmcmp/xlmass/*/lib64"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

                for dir in none $try_libdirs
                do
                        unset ac_cv_search_vexp # clear cached value
                        if test "$dir" = "none"
                        then
                                try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi
                        FFLAGS="$test_fflags"
                        LDFLAGS="$test_ldflags $try_loption"
                        LIBS=""
                        AC_SEARCH_LIBS(vexp, massvp4_64, , , -lmass_64)
                        if test "$ac_cv_search_vexp" = "-lmassvp4_64"
                        then mass_libs="$try_loption $ac_cv_search_vexp -lmass_64"
                        fi
                        if test "$ac_cv_search_vexp" != "no" ; then break ; fi
                done
                ;;

        esac
fi

if test "$mass_libs" != ""; then 
   try_dflags="$try_dflags -D__MASS"
   if test "$arch" = "ppc64-bg"; then
   # BlueGene wants this when mass libs are loaded, SP6 doesn't want this!
     ldflags="$ldflags -Wl,--allow-multiple-definition"
   fi
   if test "$arch" = "ppc64-bgq"; then
   # BlueGene wants this when mass libs are loaded, SP6 doesn't want this!
     ldflags="$ldflags -Wl,--allow-multiple-definition"
   fi
fi

echo setting MASS_LIBS... $mass_libs

# check for mpi
# some architectures require to link mpi libraries explicitly
F77=$mpif90 # use parallel compiler
if test "$mpi_libs" = ""
then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        if test "$use_parallel" -ne 0
        then
                if test "$have_mpi" -eq 0
                        # check for mpi
                then
                        unset ac_cv_search_mpi_init # clear cached value
                        LDFLAGS="$test_ldflags"
                        LIBS="$mpi_libs"
                        AC_SEARCH_LIBS(mpi_init, mpi,
                                       have_mpi=1 mpi_libs="$LIBS")
                fi
        fi
fi
echo setting MPI_LIBS... $mpi_libs


# final check on availability of parallel environment
parallel=0
for dummy in x # to allow simple 'break' 
do 
    test "$use_parallel" -eq 0 && break

    F77=$mpif90
    LIBS="$mpi_libs"
    AC_SEARCH_LIBS(mpi_init, "", parallel=1
                   try_dflags="$try_dflags -D__MPI -D__PARA")
# look for scalapack if required
    test "$parallel" -eq 0 && break
    test "$with_scalapack" -eq 0 && break
    if test "$scalapack_libs" = "" ; then
# no additional libraries needed
       AC_SEARCH_LIBS(pdgemr2d, "" , have_scalapack=1
                   try_dflags="$try_dflags -D__SCALAPACK")
       test "$have_scalapack" -eq 1 && break

# Intel MKL blacs&scalapack - Norbert Nemec 2010/08/20
# 
# Carlo Cavazzoni writes
#      If you want to use the version of scalapack distributed with MKL
#      you should be careful about:
#      1) the size of the integer variable when calling MKL subroutine,
#      QE (up to now) uses integer*4, so you should link version with
#       the suffix: lp64 and NOT ilp64 .
#      2) you should link the appropriate blacs for your MPI library:
#      INTEL MPI   -> libmkl_blacs_intelmpi_lp64.a
#      OpenMPI -> libmkl_blacs_openmpi_lp64.a
#      SGI MPI -> libmkl_blacs_sgimpt_lp64.a
#
#      if you are in doubt or if you have another version of MPI
#      you can always build BLACS by yourself.
#
if test "$have_mkl" -eq 1
   then
      unset ac_cv_search_pdgemr2d # clear cached value
      LIBS="-lmkl_blacs_lp64 $mpi_libs $blas_libs"
      if test $with_scalapack -eq 1; then
         scalapack_libs=-lmkl_blacs_openmpi_lp64
      else
         scalapack_libs=-lmkl_blacs_intelmpi_lp64
      fi
      AC_SEARCH_LIBS(pdgemr2d, "mkl_scalapack_lp64" , have_scalapack=1
                     try_dflags="$try_dflags -D__SCALAPACK"
                     scalapack_libs="-lmkl_scalapack_lp64 $scalapack_libs" )
      test "$have_scalapack" -eq 1 && break
fi
#
# sci libraries (e.g. cray xt)
       unset ac_cv_search_pdgemr2d # clear cached value
       scalapack_libs="-lsci"
       LIBS="$mpi_libs $scalapack_libs"
       AC_SEARCH_LIBS(pdgemr2d, "" , have_scalapack=1
                      try_dflags="$try_dflags -D__SCALAPACK")
       test "$have_scalapack" -eq 1 && break
# scalapack (including blacs), no -L options
       unset ac_cv_search_pdgemr2d # clear cached value
       scalapack_libs="-lscalapack"
       LIBS="$mpi_libs $scalapack_libs"
       LDFLAGS=""
       AC_SEARCH_LIBS(pdgemr2d, "" , have_scalapack=1
                   try_dflags="$try_dflags -D__SCALAPACK")
       test "$have_scalapack" -eq 1 && break
# scalapack + blacs, no -L options
       unset ac_cv_search_pdgemr2d # clear cached value
       blacs_libs="-lblacs -lblacsF77init -lblacs"
       scalapack_libs="-lscalapack $blacs_libs"
       LIBS="$mpi_libs $scalapack_libs"
       LDFLAGS=""
       AC_SEARCH_LIBS(pdgemr2d, "" , have_scalapack=1
                   try_dflags="$try_dflags -D__SCALAPACK")
       test "$have_scalapack" -eq 1 && break
# scalapack + blacs with -L options
       unset ac_cv_search_pdgemr2d # clear cached value
       if test "$scalapack_dir" = ""; then scalapack_dir="/bgsys/local/scalapack/lib"; fi
       if test "$blacs_dir" = ""; then blacs_dir="/bgsys/local/blacs/lib"; fi
       blacs_libs="-L$blacs_dir -lblacs -lblacsF77init -lblacs"
       scalapack_libs="-L$scalapack_dir -lscalapack $blacs_libs"
       LIBS="$mpi_libs $scalapack_libs"
       LDFLAGS=""
       AC_SEARCH_LIBS(pdgemr2d, "" , have_scalapack=1
                   try_dflags="$try_dflags -D__SCALAPACK")
    else
        # scalapack provided in SCALAPACK_LIBS - not checked!
        have_scalapack=1
        try_dflags="$try_dflags -D__SCALAPACK"
    fi   
done

# ELPA iff SCALAPACK (statically linked)
elpa_libs_switch="disabled"
if test "$with_elpa" -eq 2 && test "$have_scalapack" -eq 1; then
    if test "$use_openmp" -eq 1 ; then
        scalapack_libs="$with_elpa_path/lib/libelpa_mt.a $scalapack_libs"
        try_iflags="$try_iflags -I$with_elpa_path/include/elpa/modules "
        try_dflags="$try_dflags -D__ELPA"
    else
        scalapack_libs="$with_elpa_path/libelpa.a $scalapack_libs"
        try_iflags="$try_iflags -I$with_elpa_path/include/elpa/modules "
        try_dflags="$try_dflags -D__ELPA"
    fi
fi

if test "$with_elpa" -eq 1 && test "$have_scalapack" -eq 1; then
    scalapack_libs="\$(TOPDIR)/ELPA/libelpa.a $scalapack_libs"
    try_dflags="$try_dflags -D__ELPA"
    elpa_libs_switch="enabled"    
fi
    
echo setting SCALAPACK_LIBS... $scalapack_libs

# preprocessing flag for openmp (experimental)
if test "$use_openmp" -eq 1 ; then try_dflags="$try_dflags -D__OPENMP" ; fi

# preprocessing flag for signal trapping (experimental)
if test "$use_signals" -eq 1 ; then try_dflags="$try_dflags -D__TRAP_SIGUSR1" ; fi

#OBM - This is yet another variable for internal copy of BLAS/LAPACK
if test "$extlib_flags" = "" ; then 
       case "$arch" in
        ppc64 )
		extlib_flags="-q64 -qthreaded"
		echo setting extlib_flags ... $extlib_flags
        ;;
       esac
fi
##


if test "$dflags" = "" ; then dflags="$try_dflags" ; fi
echo setting DFLAGS... $dflags

if test "$iflags" = "" ; then iflags="$try_iflags" ; fi
echo setting IFLAGS... $iflags

# xlf compilers (AIX and powerpc) want comma-separated -D directives
if test "$xlf_flags" -ne 0
then
        fdflags="`echo $dflags | sed 's/  */,/g'`"
else
        fdflags="\$(DFLAGS)"
fi
echo setting FDFLAGS... $fdflags

if test "$ranlib" != "echo"
then
   AC_CHECK_PROG(ranlib,ranlib,ranlib,echo)
fi
# MacOs workaround for ranlib
if test "$arch" = "mac686"; then
   if test "$ranlib" = "ranlib"; then
      ranlib="ranlib -c"
   fi
fi

echo setting RANLIB... $ranlib

# check if wget or curl work  
AC_CHECK_PROG(wget, wget, wget -O)
if test "$wget" = ""; then
   AC_CHECK_PROG(wget, curl, curl -o)
fi
echo setting WGET... $wget

# configure output messages
blas_line="BLAS_LIBS=$blas_libs"
if test "$have_scalapack" -eq 1; then
   scalapack_line="SCALAPACK_LIBS=$scalapack_libs"
else
   scalapack_libs=""
   scalapack_line="@delete@"
fi
lapack_line="LAPACK_LIBS=$lapack_libs"
fft_line="FFT_LIBS=$fft_libs"
if test "$mpi_libs" != "" ; then
   mpi_line="MPI_LIBS=$mpi_libs"
else
   mpi_line="@delete@"
fi
if test "$mass_libs" != "" ; then
   mass_line="MASS_LIBS=$mass_libs"
else
   mass_line="@delete@"
fi

if test "$use_parallel" -ne 0
then
        if test "$parallel" -ne 0
        then
                parallel_report="Parallel environment detected successfully.\\
Configured for compilation of parallel executables."
        else
                parallel_report="Parallel environment not detected \
\(is this a parallel machine?\).\\
Configured for compilation of serial executables."
        fi
else
        parallel_report="Configured for compilation of serial executables."
fi

# export settings to generated files
AC_SUBST(cc)
AC_SUBST(cflags)
AC_SUBST(dflags)
AC_SUBST(fdflags)
AC_SUBST(cpp)
AC_SUBST(cppflags)
AC_SUBST(f90)
AC_SUBST(mpif90)
AC_SUBST(f90flags)
AC_SUBST(f77)
AC_SUBST(fflags)
AC_SUBST(fflags_noopt)
AC_SUBST(fflags_nomain)
AC_SUBST(pre_fdflags)
AC_SUBST(imod)
AC_SUBST(iflags)
AC_SUBST(ld)
AC_SUBST(ldflags)
AC_SUBST(ld_libs)
AC_SUBST(blas_libs)
AC_SUBST(blas_libs_switch)
AC_SUBST(lapack_libs)
AC_SUBST(lapack_libs_switch)
AC_SUBST(fft_libs)
AC_SUBST(mpi_libs)
AC_SUBST(mass_libs)
AC_SUBST(elpa_libs_switch)
AC_SUBST(scalapack_libs)
AC_SUBST(ar)
AC_SUBST(arflags)
AC_SUBST(ranlib)
AC_SUBST(f90rule)
AC_SUBST(blas_line)
AC_SUBST(lapack_line)
AC_SUBST(scalapack_line)
AC_SUBST(fft_line)
AC_SUBST(mpi_line)
AC_SUBST(mass_line)
AC_SUBST(parallel_report)
AC_SUBST(wget)
AC_SUBST(topdir)
AC_SUBST(extlib_flags)

AC_CONFIG_HEADERS([include/c_defs.h:include/c_defs.h.in])
AC_CONFIG_FILES(make.sys)
AC_CONFIG_FILES(configure.msg)

AC_CONFIG_FILES(install/make_wannier90.sys)

AC_CONFIG_FILES(install/make_blas.inc)
AC_CONFIG_FILES(install/make_lapack.inc)


AC_OUTPUT

# final messages
sed '/@delete@/d' configure.msg
echo configure: success
