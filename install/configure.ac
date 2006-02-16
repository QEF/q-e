AC_INIT(ESPRESSO, 3.0, , espresso)

# -----------------------------
# application-specific settings

# default choices
# set each variable to 0 (false) or 1 (true)
use_parallel=1
use_shared=1

# -----------------------------

# store variables from the environment (may be set or not)
topdir=$TOPDIR
arch=$ARCH
env_cc=$CC
mpicc=$MPICC
cpp=$CPP
cflags=$CFLAGS
cppflags=$CPPFLAGS
dflags=$DFLAGS
f77=$F77
f90=$F90
mpif77=$MPIF77
mpif90=$MPIF90
fflags=$FFLAGS
f90flags=$F90FLAGS
ld=$LD
ldflags=$LDFLAGS
blas_libs=$BLAS_LIBS
lapack_libs=$LAPACK_LIBS
fft_libs=$FFT_LIBS
mpi_libs=$MPI_LIBS
mass_libs=$MASS_LIBS
libdirs=$LIBDIRS
ar=$AR
arflags=$ARFLAGS
mylib=$MYLIB
includefftw=$INCLUDEFFTW

# command-line arguments
for arg in $*
do
	if   test "$arg" = "--enable-parallel"  ; then use_parallel=1
	elif test "$arg" = "--disable-parallel" ; then use_parallel=0
	elif test "$arg" = "--enable-shared"    ; then use_shared=1
	elif test "$arg" = "--disable-shared"   ; then use_shared=0
	fi
done

# configure for current directory by default
if test "$topdir" = "" ; then topdir="`pwd`" ; fi

# check system type (no cross-compilation for now)
AC_CANONICAL_BUILD

# many HPC systems are configured so that running parallel programs
# interactively is disabled: on those systems, AC_PROG_F77 and AC_PROG_CC
# would fail because they can't run the compiled executables
# to work around that, let's pretend we are cross-compiling even if we aren't
# !!! this relies on undocumented Autoconf behavior !!! 
if test "$host" = "" ; then host=$build ; fi
cross_compiling=yes

# identify architecture
if test "$arch" = ""
then
	case $host in
		ia64-*-linux-gnu )      arch=ia64   ;;
		x86_64-*-linux-gnu )    arch=amd64  ;;
		*-pc-linux-gnu )        arch=ia32   ;;
		*-ibm-aix* )            arch=aix    ;;
		mips-sgi-* )            arch=mips   ;;
		alphaev*-dec-osf* )     arch=alpha  ;;
		alphaev*-*-linux* )     arch=alinux ;;
		sparc-sun-* )           arch=sparc  ;;
		*cray-unicosmp* )       arch=crayx1 ;;
		powerpc-apple-darwin* ) arch=mac    ;;
		*-pc-cygwin )           arch=cygwin ;;
	esac
fi
echo checking architecture... $arch

# candidate compilers, archiver and flags based on architecture
case $arch in
ia64 | amd64 )
	try_mpif90="mpif90"
	try_f90="ifort ifc efc pgf90 g95 f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__LINUX64"
	;;
ia32 )
	try_mpif90="mpif90"
	try_f90="ifort ifc pgf90 g95 f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__LINUX"
	;;
aix )
	try_mpif90="mpxlf90_r mpxlf90"
	try_f90="xlf90_r xlf90 f90"
	try_ar="ar"
	try_arflags="-X64 ruv"
	try_dflags="-D__AIX -D__XLF"
	;;
mips )
	try_f90="f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__COMPLIB -D__SGI64 -D__ORIGIN"
	;;
alpha )
	try_f90="f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__ALPHA"
	;;
alinux )
	try_mpif90="mpif90"
	try_f90="fort g95 f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__ALPHA -D__LINUX64"
	;;
sparc )
	try_mpif90="mpf90"
	try_f90="f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__SUN"
	;;
crayx1 )
	try_f90="ftn"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__X1"
	;;
mac )
	try_mpif90="mpxlf90 mpf90"
	try_f90="xlf90 f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__MAC -D__XLF"
	ranlib="ranlib"
	;;
cygwin )
	try_mpif90="mpif90"
	try_f90="g95 f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__LINUX"
	;;
* )
	AC_MSG_WARN(unsupported architecture)
	try_f90="f90"
	try_ar="ar"
	try_arflags="ruv"
	;;
esac

if test "$use_parallel" -ne 0
then
	# clear cached values
	unset F77 ac_cv_prog_ac_ct_F77 ac_cv_f77_compiler_gnu ac_cv_prog_f77_g

	# check parallel Fortran 90 compiler
	if test "$mpif90" = "" -a "$f90" != "" ; then mpif90=$f90; fi
	if test "$mpif90" = "" ; then mpif90="$try_mpif90 $try_f90" ; fi
	AC_PROG_F77($mpif90)
	mpif90=$F77

	case "$arch" in
	ia32 | ia64 | amd64 )
		# check which compiler does mpif90 wrap
		echo $ECHO_N "checking version of $mpif90... $ECHO_C"
		ifort_version=`$mpif90 -V 2>&1 | grep "Intel(R)"`
		pgf_version=`$mpif90 -V 2>&1 | grep "^pgf"`
		g95_version=`$mpif90 -v 2>&1 | grep "g95"`
		if test "$ifort_version" != ""
		then
			version=`$mpif90 -V 2>&1 | grep Version |
				 sed 's/.*Version//' | awk '{print $1}'`
			if test `echo $version | sed 's/\..*//'` -lt 8
			then
				try_f90="ifc efc"
				echo "${ECHO_T}ifc $version"
			else
				try_f90="ifort efc"
				echo "${ECHO_T}ifort $version"
			fi
		elif test "$pgf_version" != ""
		then
			version=`echo $pgf_version | awk '{print $2}'`
			try_f90="pgf90"
			echo "${ECHO_T}pgf90 $version"
		elif test "$g95_version" != ""
		then
			version=`echo $g95_version | awk '{print $3}'`
			try_f90="g95"
			echo "${ECHO_T}g95 $version"
		else
			echo "${ECHO_T}unknown"
		fi
		;;
	esac

	case "$arch:$mpif90" in
	*:mpif90 )
		try_mpicc="mpicc"
		try_mpif77="mpif77"
		;;
	*:mpf90 )
		try_mpicc="mpcc"
		try_mpif77="mpf77"
		;;
	aix:mpxlf*_r )
		try_mpicc="mpcc_r"
		try_mpif77="mpxlf_r"
		;;
	aix:mpxlf* )
		try_mpicc="mpcc"
		try_mpif77="mpxlf"
		;;
	mac:mpxlf* )
		try_mpicc="mpcc"
		try_mpif77="mpxlf mpf77"
		;;
	esac
fi

# clear cached values
unset F77 ac_cv_prog_ac_ct_F77 ac_cv_f77_compiler_gnu ac_cv_prog_f77_g

# check serial Fortran 90 compiler
if test "$f90" = "" ; then f90="$try_f90" ; fi
AC_PROG_F77($f90)
f90=$F77
if test "$use_parallel" -eq 0 ; then mpif90=$f90 ; fi

echo setting F90... $f90
echo setting MPIF90... $mpif90

case "$arch:$f90" in
*:ifort* )
	try_cc="icc ecc gcc cc"
	try_f77="$f90"
	;;
*:ifc* )
	try_cc="icc ecc gcc cc"
	try_f77="$f90"
	;;
*:pgf90 )
	try_cc="pgcc gcc cc"
	try_f77="pgf77 $f90"
	;;
*:g95 )
	try_cc="gcc cc"
	try_f77="$f90"
	;;
ia32:* | ia64:* | amd64:* )
	try_cc="gcc cc"
	try_f77="f77 $f90"
	;;
aix:*xlf*_r )
	try_cc="xlc_r cc gcc"
	try_f77="xlf_r f77 $f90"
	;;
aix:*xlf* )
	try_cc="xlc cc gcc"
	try_f77="xlf f77 $f90"
	;;
mips:* )
	try_cc="cc gcc"
	try_f77="f77 $f90"
	;;
alpha:* )
	try_cc="cc gcc"
	try_f77="f77 $f90"
	;;
alinux:* )
	try_cc="ccc gcc"
	try_f77="fort g77 f77 $f90"
	;;
sparc:* )
	try_cc="cc gcc"
	try_f77="f77 $f90"
	;;
crayx1:ftn )
	try_cc="cc"
	try_f77="$f90"
	;;
mac:xlf90 )
	try_cc="gcc"
	try_f77="xlf f77 $f90"
	;;
mac:* )
	try_cc="gcc"
	try_f77="f77 $f90"
	;;
cygwin*:g95 )
	try_cc="gcc cc"
	try_f77="g95 $f90"
	;;
* )
	# unknown, try these
	try_cc="cc gcc"
	try_f77="f77 $f90"
	;;
esac

if test "$use_parallel" -ne 0
then
	# clear cached values
	unset CC ac_cv_prog_ac_ct_CC ac_cv_c_compiler_gnu ac_cv_prog_cc_g

	# check parallel C compiler
	if test "$mpicc" = "" ; then mpicc="$try_mpicc $try_cc" ; fi
	AC_PROG_CC($mpicc)
	mpicc=$CC

	# clear cached values
	unset F77 ac_cv_prog_ac_ct_F77 ac_cv_f77_compiler_gnu ac_cv_prog_f77_g

	# check parallel Fortran 77 compiler
	if test "$mpif77" = "" ; then mpif77="$try_mpif77 $try_f77" ; fi
	AC_PROG_F77($mpif77)
	mpif77=$F77
fi

# clear cached values
unset CC ac_cv_prog_ac_ct_CC ac_cv_c_compiler_gnu ac_cv_prog_cc_g

# check serial C compiler
if test "$env_cc" = "" ; then cc="$try_cc" ; else cc="$env_cc"; fi
AC_PROG_CC($cc)
cc=$CC
if test "$use_parallel" -eq 0 ; then mpicc=$cc ; fi

echo setting CC... $cc
echo setting MPICC... $mpicc

# clear cached values
unset F77 ac_cv_prog_ac_ct_F77 ac_cv_f77_compiler_gnu ac_cv_prog_f77_g

# check serial Fortran 77 compiler
if test "$f77" = "" ; then f77="$try_f77" ; fi
AC_PROG_F77($f77)
f77=$F77
if test "$use_parallel" -eq 0 ; then mpif77=$f77 ; fi

echo setting F77... $f77
echo setting MPIF77... $mpif77

# check whether the C and Fortran compilers are compatible
case "$arch:$cc:$f90" in
	# list supported combinations here
ia32:icc:ifort* | ia32:icc:ifc* )   ;;
ia32:gcc:ifort* | ia32:gcc:ifc* )   ;;
ia32:pgcc:pgf90 )                   ;;
ia32:gcc:g95 )                      ;;
ia64:icc:ifort* | ia64:icc:ifc* )   ;;
ia64:gcc:ifort* | ia64:gcc:ifc* )   ;;
ia64:pgcc:pgf90 )                   ;;
ia64:gcc:g95 )                      ;;
amd64:icc:ifort* | amd64:icc:ifc* ) ;;
amd64:gcc:ifort* | amd64:gcc:ifc* ) ;;
amd64:pgcc:pgf90 )                  ;;
amd64:gcc:g95 )                     ;;
aix:mpcc*:mpxlf* )                  ;;
aix:xlc*:*xlf* | aix:cc:*xlf* )     ;;
mips:cc:f90 )                       ;;
alinux:ccc:fort | alinux:gcc:fort ) ;;
alpha:cc:f90 )                      ;;
mac:mpcc:mpf90 )                    ;;
mac:gcc:xlf90 )                     ;;
crayx1:cc:ftn )                     ;;
sparc:cc:f90 )                      ;;
cygwin:gcc:g95 )                    ;;

* )
	AC_MSG_WARN(unsupported C/Fortran compilers combination)
	;;
esac

# check Fortran compiler flags
have_cpp=0
xlf_flags=0

case "$arch:$f90" in
ia64:ifort* | amd64:ifort* )
	try_fflags="-O2 -assume byterecl"
	try_f90flags="\$(FFLAGS) -nomodule"
	try_fflags_noopt="-O0 -assume byterecl"
	try_f90flags_noopt="\$(FFLAGS_NOOPT) -nomodule"
	try_ldflags=""
	try_ldflags_static="-static"
	try_dflags="$try_dflags -D__INTEL"
	pre_fdflags="-fpp "
	have_cpp=1
	;;
ia64:ifc* )
	try_fflags="-Vaxlib -O2"
	try_f90flags="\$(FFLAGS) -nomodule"
	try_fflags_noopt="-O0"
	try_f90flags_noopt="\$(FFLAGS_NOOPT) -nomodule"
	try_ldflags="-Vaxlib"
	try_ldflags_static="-static"
	try_dflags="$try_dflags -D__INTEL"
	pre_fdflags="-fpp "
	have_cpp=1
	;;
ia64:pgf* | amd64:pgf* )
	try_fflags="-fast -r8"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0"
	try_f90flags_noopt="\$(FFLAGS_NOOPT)"
	try_ldflags=""
	try_ldflags_static="-Bstatic"
	try_dflags="$try_dflags -D__PGI"
	;;
ia32:ifort* )
        try_fflags="-O2 -tpp6 -assume byterecl"
	try_f90flags="\$(FFLAGS) -nomodule"
	try_fflags_noopt="-O0 -assume byterecl"
	try_f90flags_noopt="\$(FFLAGS_NOOPT) -nomodule"
	try_ldflags=""
	try_ldflags_static="-static"
	try_dflags="$try_dflags -D__INTEL"
	pre_fdflags="-fpp "
	have_cpp=1
	;;
ia32:ifc* )
        try_fflags="-Vaxlib -O2 -tpp6"
	try_f90flags="\$(FFLAGS) -nomodule"
	try_fflags_noopt="-O0"
	try_f90flags_noopt="\$(FFLAGS_NOOPT) -nomodule"
	try_ldflags="-Vaxlib"
	try_ldflags_static="-static"
	try_dflags="$try_dflags -D__INTEL"
	pre_fdflags="-fpp "
	have_cpp=1

	if test `echo $version | sed 's/\..*//'` -lt 7
	then
		# old versions of ifc require this stuff
		echo $ECHO_N "setting up ifc environment... $ECHO_C"
		try_f90flags="\$(FFLAGS) -cl,./intel.pcl"
		try_f90flags_noopt="\$(FFLAGS_NOOPT) -cl,./intel.pcl"
		pcl_ph="D3 Raman"
		pcl_pw="PH PP Gamma PWCOND pwtools $pcl_ph"
		pcl_modules="PW CPV flib upftools atomic $pcl_pw"
		pcl_dot=". Modules $pcl_modules"
		for dir in $pcl_dot
		do
			echo work.pc > $topdir/$dir/intel.pcl
		done
		for dir in $pcl_modules
		do
			echo ../Modules/work.pc >> $topdir/$dir/intel.pcl
		done
		for dir in $pcl_pw
		do
			echo ../PW/work.pc >> $topdir/$dir/intel.pcl
		done
		for dir in $pcl_ph
		do
			echo ../PH/work.pc >> $topdir/$dir/intel.pcl
		done
		echo "${ECHO_T}done"
	fi
	;;
ia32:pgf* )
	try_fflags="-fast -r8"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0"
	try_f90flags_noopt="\$(FFLAGS_NOOPT)"
	try_ldflags=""
	try_ldflags_static="-Bstatic"
	try_dflags="$try_dflags -D__PGI"
	;;
*:g95 )
	try_fflags="-O3 -cpp"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0 -cpp"
	try_f90flags_noopt="\$(FFLAGS_NOOPT)"
	try_ldflags=""
	#try_ldflags_static="-static"
	try_dflags="$try_dflags -D__G95"
	have_cpp=1
	;;
cygwin:g95 )
        try_fflags="-O3 -cpp"
        try_f90flags="\$(FFLAGS)"
        try_fflags_noopt="-O0 -cpp"
        try_f90flags_noopt="\$(FFLAGS_NOOPT)"
        try_ldflags=""
        #try_ldflags_static="-static"
        try_dflags="$try_dflags -D__G95"
        have_cpp=1
        ;;
aix:*xlf* )
	try_fflags="-q64 -qalias=noaryovrlp -O3 -qstrict \
-qarch=auto -qtune=auto -qdpc -Q -qalias=nointptr"
	try_f90flags="\$(FFLAGS) -qsuffix=cpp=f90 -qfree=f90"
	try_fflags_noopt="-q64 -O0"
	try_f90flags_noopt="\$(FFLAGS_NOOPT) -qsuffix=cpp=f90 -qfree=f90"
	try_ldflags="-q64"
	# try_ldflags_static="-bstatic"
	pre_fdflags="-WF,"
	xlf_flags=1
	have_cpp=1
	;;
mips:f90 )
	try_fflags="-mips4 -64 -O2 -r10000 -r8"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-mips4 -64 -O0"
	try_f90flags_noopt="\$(FFLAGS_NOOPT)"
	try_ldflags="-mips4 -64"
	pre_fdflags="-ftpp -macro_expand "
	have_cpp=1
	;;
alinux:fort )
	have_cpp=1
	pre_fdflags="-cpp "
	try_fflags="-O -r8 -align dcommons -align records"
	try_f90flags="\$(FFLAGS) -free"
	try_fflags_noopt="-O0 -r8 -align dcommons -align records"
	try_f90flags_noopt="\$(FFLAGS_NOOPT) -free"
        try_ldflags_static="-non_shared"
	;;
alpha:f90 )
	have_cpp=1
	pre_fdflags="-cpp "
	try_fflags="-O -real_size 64 -align dcommons -align records -assume buffered_io"
	try_f90flags="\$(FFLAGS) -free"
	try_fflags_noopt="-O0 -real_size 64 -align dcommons -align records -assume buffered_io"
	try_f90flags_noopt="\$(FFLAGS_NOOPT) -free"
	;;
sparc:mpf90 | sparc:f90 )
	try_fflags="-fast -O1 -nodepend -xvector=no -xchip=ultra3 \
-xarch=v8plusb -xlic_lib=sunperf"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0 -xlic_lib=sunperf"
	try_f90flags_noopt="\$(FFLAGS_NOOPT)"
	try_ldflags=""
	imod="-M"
	;;
crayx1:ftn )
        try_fflags="-s default64 -dp -rma -e0 \
-O scalar3,stream3,vector3,nointerchange -Z -O inline5"
        try_f90flags="\$(FFLAGS)"
        try_fflags_noopt="-s default64 -dp -rma -e0 -O0 -Z"
        try_f90flags_noopt="\$(FFLAGS_NOOPT)"
        try_ldflags="-s default64 -f nan64"
        try_dflags="$try_dflags -D__X1_COA"
        pre_fdflags="-e Z -F"
        have_cpp=1
	;;
mac:* )
	try_fflags="-O4 -qarch=auto -qtune=auto -qsuffix=cpp=f90 \
-qdpc -qalias=nointptr"
	try_f90flags="\$(FFLAGS) -qfree=f90"
	try_fflags_noopt="-O0"
	try_f90flags_noopt="\$(FFLAGS_NOOPT) -qfree=f90"
	try_ldflags=""
	pre_fdflags="-WF,"
	xlf_flags=1
	have_cpp=1
	;;
* )
	# unknown, try these
	try_fflags="-O"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0"
	try_f90flags_noopt="$\(FFLAGS_NOOPT)"
	try_ldflags=""
	;;
esac
if test "$use_shared" -eq 0 ; then
try_ldflags="$try_ldflags $try_ldflags_static" ; fi

if test "$fflags" = ""   ; then fflags=$try_fflags     ; fi
if test "$f90flags" = "" ; then f90flags=$try_f90flags ; fi
if test "$fflags_noopt" = ""   ; then fflags_noopt=$try_fflags_noopt     ; fi
if test "$f90flags_noopt" = "" ; then f90flags_noopt=$try_f90flags_noopt ; fi
echo setting FFLAGS... $fflags
echo setting F90FLAGS... $f90flags
echo setting FFLAGS_NOOPT... $fflags_noopt
echo setting F90FLAGS_NOOPT... $f90flags_noopt
if test "$imod" = "" ; then imod="-I" ; fi

case "$arch:$cc" in
*:icc )
	try_cflags="-O3"
	c_ldflags=""
	;;
*:pgcc )
	try_cflags="-fast"
	c_ldflags=""
	try_cpp="cpp"
	;;
cygwin:gcc )
        try_cflags="-O3 -fomit-frame-pointer"
        c_ldflags=""
        try_cpp="cpp"
        ;;
aix:mpcc* | aix:xlc* | aix:cc )
	try_cflags="-q64 -O2"
	c_ldflags="-q64"
	;;
mips:cc )
	try_cflags="-mips4 -64 -O2 -r10000"
	;;
alpha:cc )
	try_cflags="-O"
	;;
sparc:mpcc | sparc:cc )
	try_cflags="-fast -dalign -xchip=ultra3 -xarch=v8plusb \
-xlic_lib=sunperf"
	try_cpp="fpp"
	;;
crayx1:cc )
        try_cflags=""
        c_ldflags=""
        ;;
mac:gcc )
	try_cflags="-O3 -fomit-frame-pointer -I/usr/include/malloc"
	try_cpp="cpp"
	;;
mac:* )
	try_cflags="-O4"
	;;
*:gcc )
	try_cflags="-O3 -fomit-frame-pointer"
	c_ldflags=""
	try_cpp="cpp"
	;;
* )
	try_cflags="-O"
	;;
esac
if test "$cflags" = "" ; then cflags=$try_cflags ; fi
echo setting CFLAGS... $cflags

# preprocessor
AC_PROG_CPP
if test "$cpp" = "" ; then cpp=$try_cpp; fi
if test "$cpp" = "" ; then cpp=$CPP;     fi
echo setting CPP... $cpp

echo $ECHO_N "setting CPPFLAGS... $ECHO_C"
case $cpp in
	cpp) try_cppflags="-P -traditional" ;;
	fpp) try_cppflags="-P"              ;;
	*)   try_cppflags=""                ;;
esac
if test "$cppflags" = "" ; then cppflags=$try_cppflags ; fi
echo "${ECHO_T}$cppflags"

# linker and archiver
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
	f90rule="\$(CPP) \$(CPPFLAGS) \$< -o \$*.F90 \\
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

AC_LANG_PUSH(Fortran 77)
F77=$f90 # use Fortran 90 actually (autoconf only knows Fortran 77!)

# search for libraries
have_blas=0
have_mkl=0
have_lapack=0
have_essl=0
have_fft=0
have_fftw=0
have_fftw_h=0
have_mpi=0

# check for blas
# supported vendor replacements:
#   mkl on linux
#   essl on aix
#   complib.sgimath on mips
#   cxml on alpha
#   SUNperf on sparc
# atlas is used over blas if available
# internal version is used if none is found
if test "$blas_libs" = ""
then
	# check directories in LD_LIBRARY_PATH too
	# (maybe they are already searched by default, but I'm not sure)
	ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

	# build list of directories to search
	case "$arch" in
	mac )    try_libdirs="/usr/local/lib /sw/lib" ;;
	*)       try_libdirs="/usr/local/lib"         ;;
	esac

	# search
	case "$arch:$f90" in

	amd64:pgf* )
		# check for acml
		unset ac_cv_search_dgemm # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS=""
		AC_SEARCH_LIBS(dgemm, acml, have_blas=1 blas_libs="$LIBS")
		;;

	*:pgf* )
		# check for PGI blas
		unset ac_cv_search_dgemm # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS=""
		AC_SEARCH_LIBS(dgemm, blas, have_blas=1 blas_libs="$LIBS")
		;;

	ia64:* )
		# check for mkl (in several directories)
		try_libdirs="/opt/intel/mkl/8.0/lib/64
			     /opt/intel/mkl721/lib/64
			     /opt/intel/mkl72/lib/64
			     /opt/intel/mkl701/lib/64
			     /opt/intel/mkl70/lib/64
			     /opt/intel/mkl/mkl61/lib/64
			     /opt/intel/mkl/lib/64
			     /opt/intel/mkl61/lib/64"
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
			LIBS="-lguide -lpthread"
			AC_SEARCH_LIBS(dgemm, mkl_ipf, have_blas=1 have_mkl=1
				       blas_libs="$try_loption $LIBS")
			if test "$ac_cv_search_dgemm" != "no"
			then break ; fi
		done
		;;

	amd64:* )
		try_libdirs="/opt/intel/mkl/8.0/lib/em64t
			     /opt/intel/mkl721/lib/em64t"
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
			LIBS="-lguide -lpthread"
			AC_SEARCH_LIBS(dgemm, mkl_em64t, have_blas=1 have_mkl=1
				       blas_libs="$try_loption $LIBS")
			if test "$ac_cv_search_dgemm" != "no"
			then break ; fi
		done
		;;

	ia32:* )
		# check for mkl (in several directories)
		try_libdirs="/opt/intel/mkl/8.0/lib/32
			     /opt/intel/mkl721/lib/32
			     /opt/intel/mkl72/lib/32
			     /opt/intel/mkl701/lib/32
			     /opt/intel/mkl70/lib/32
			     /opt/intel/mkl/mkl61/lib/32
			     /opt/intel/mkl/lib/32
			     /opt/intel/mkl61/lib/32
			     /cineca/prod/intel/lib"
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
			LIBS="-lguide -lpthread"
			AC_SEARCH_LIBS(dgemm, mkl_ia32, have_blas=1 have_mkl=1
				       blas_libs="$try_loption $LIBS")
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
		AC_SEARCH_LIBS(dgemm, essl, have_essl=1 have_blas=1
			       blas_libs="$LIBS")
		;;

	mips:* )
		# check for complib.sgimath
		unset ac_cv_search_dgemm # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS=""
		AC_SEARCH_LIBS(dgemm, complib.sgimath, have_blas=1
			       blas_libs="$LIBS")
		;;

	alpha:* | alinux:* )
		# check for cxml
		unset ac_cv_search_dgemm # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS=""
		AC_SEARCH_LIBS(dgemm, cxml, have_blas=1 have_lapack=1
                                blas_libs="$LIBS")
		;;

	sparc:* )
		# check for SUNperf library
		unset ac_cv_search_zggev # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS=""
		AC_SEARCH_LIBS(dgemm, sunperf, have_blas=1
			       blas_libs="-xlic_lib=sunperf $LIBS")
		;;

	esac

	if test "$have_blas" -eq 0
	then
		# check for atlas (in several directories)
		case "$f90" in
		ifort* | ifc* )
			try_libdirs="/usr/local/lib /cineca/prod/intel/lib" ;;
		pgf* )	try_libdirs="/usr/local/lib /cineca/prod/pgi/lib"   ;;
		g95 )	try_libdirs="/usr/local/lib /cineca/prod/gnu/lib"   ;;
		* )	try_libdirs="/usr/local/lib"                        ;;
		esac
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
			LIBS="-latlas -lg2c"
			AC_SEARCH_LIBS(dgemm, f77blas, have_blas=1
				       blas_libs="$try_loption $LIBS")
			if test "$ac_cv_search_dgemm" != "no"
			then break ; fi
		done
	fi

	if test "$have_blas" -eq 0
	then
		# check for blas (in several directories)
		case "$f90" in
		ifort* | ifc* )
			try_libdirs="/usr/local/lib /cineca/prod/intel/lib" ;;
		pgf* )	try_libdirs="/usr/local/lib /cineca/prod/pgi/lib"   ;;
		g95 )	try_libdirs="/usr/local/lib /cineca/prod/gnu/lib"   ;;
		* )	try_libdirs="/usr/local/lib"                        ;;
		esac
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
			AC_SEARCH_LIBS(dgemm, blas, have_blas=1
				       blas_libs="$try_loption $LIBS")
			if test "$ac_cv_search_dgemm" != "no"
			then break ; fi
		done
	fi
fi
echo setting BLAS_LIBS... $blas_libs

# check for lapack
# supported vendor replacements:
#   mkl on linux
#   SUNperf on sparc
# internal version is used if none is found
if test "$lapack_libs" = ""
then
	# check directories in LD_LIBRARY_PATH too
	# (maybe they are already searched by default, but I'm not sure)
	ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

	# build list of directories to search
	case "$arch" in
	mac )    try_libdirs="/usr/local/lib /sw/lib" ;;
	*)       try_libdirs="/usr/local/lib"         ;;
	esac

	# search
	case "$arch:$f90" in

	ia32:* | ia64:* | amd64:* )
		# check for mkl_lapack
		if test "$have_mkl" -ne 0
		then
			unset ac_cv_search_zggev # clear cached value
			FFLAGS="$test_fflags"
			LDFLAGS="$test_ldflags"
			LIBS="$blas_libs"
			AC_SEARCH_LIBS(zggev, mkl_lapack, have_lapack=1
				       lapack_libs="-lmkl_lapack")
		fi
		;;

	sparc:* )
		# check for SUNperf library
		unset ac_cv_search_zggev # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS="$blas_libs"
		AC_SEARCH_LIBS(zggev, sunperf, have_lapack=1
			       lapack_libs="-xlic_lib=sunperf $LIBS")
		;;

	esac

	if test "$have_lapack" -eq 0
	then
		# check for lapack (in several directories)
		case "$f90" in
		ifort* | ifc* )
			try_libdirs="/usr/local/lib /cineca/prod/intel/lib" ;;
		pgf* )	try_libdirs="/usr/local/lib /cineca/prod/pgi/lib"   ;;
		g95 )	try_libdirs="/usr/local/lib /cineca/prod/gnu/lib"   ;;
		*xlf* ) try_libdirs="/usr/local/lib /cineca/lib"            ;;
		* )	try_libdirs="/usr/local/lib"                        ;;
		esac
		try_libdirs="$libdirs $try_libdirs $ld_library_path"

		for dir in none $try_libdirs
		do
			unset ac_cv_search_zggev # clear cached value
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
			AC_SEARCH_LIBS(zggev, lapack, have_lapack=1
				       lapack_libs="$try_loption -llapack")
			if test "$ac_cv_search_zggev" != "no"
			then break ; fi
		done
	fi
	# essl must precede lapack (if present)
	if test "$have_essl" -ne 0 && test "$have_lapack" -ne 0
	then lapack_libs="-lessl $lapack_libs" ; fi
fi
echo setting LAPACK_LIBS... $lapack_libs

# check for fftw
# supported vendor replacements:
#   essl on aix
#   SUNperf on sparc
# internal version is used if none is found
if test "$fft_libs" = ""
then
	# check directories in LD_LIBRARY_PATH too
	# (maybe they are already searched by default, but I'm not sure)
	ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

	# build list of directories to search
	case "$arch" in
	mac )    try_libdirs="/usr/local/lib /sw/lib" ;;
	*)       try_libdirs="/usr/local/lib"         ;;
	esac

	# search
	case "$arch" in
	aix )
		# check for essl
		unset ac_cv_search_dcft # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS="$fft_libs"
		AC_SEARCH_LIBS(dcft, essl, have_fft=1 fft_libs="$LIBS")
		;;
	sparc )
		# check for SUNperf library
		unset ac_cv_search_zfft3i # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS="$libs"
		AC_SEARCH_LIBS(zfft3i, sunperf, have_fft=1
			       fft_libs="-xlic_lib=sunperf $LIBS")
		;;
	esac

	if test "$have_fft" -eq 0
	then
		AC_LANG_PUSH(C)

		# check for fftw (in several directories)
		case "$f90" in
		ifort* | ifc* )
			try_libdirs="/usr/local/lib /cineca/prod/intel/lib" ;;
		pgf* )	try_libdirs="/usr/local/lib /cineca/prod/pgi/lib"   ;;
		* )	try_libdirs="/usr/local/lib /cineca/prod/gnu/lib"   ;;
		esac
		try_libdirs="$libdirs $try_libdirs $ld_library_path"

		for dir in none $try_libdirs
		do
			unset ac_cv_search_fftwnd # clear cached value
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
			AC_SEARCH_LIBS(fftwnd, fftw, have_fftw=1
				       fft_libs="$try_loption $LIBS", , -lm)
			if test "$ac_cv_search_fftwnd" != "no"
			then break ; fi
		done

		AC_LANG_POP(C)
	fi
fi
echo setting FFT_LIBS... $fft_libs

# check for mass on aix
if test "$mass_libs" = ""
then
	# check directories in LD_LIBRARY_PATH too
	# (maybe they are already searched by default, but I'm not sure)
	ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

	# build list of directories to search
	case "$arch" in
	mac )    try_libdirs="/usr/local/lib /sw/lib" ;;
	*)       try_libdirs="/usr/local/lib"         ;;
	esac

	# search
	case "$arch" in
	aix )
		# check for mass (in several directories)
		try_libdirs="/usr/local/lib /cineca/lib /cineca/lib/mass"
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
				-or "$ac_cv_search_vexp" = "-lmassv"
			then mass_libs="$try_loption $ac_cv_search_vexp -lmass"
			fi
			if test "$ac_cv_search_vexp" != "no" ; then break ; fi
		done
		;;
	esac
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

	# build list of directories to search
	case "$arch" in
	mac )    try_libdirs="/usr/local/lib /sw/lib" ;;
	*)       try_libdirs="/usr/local/lib"         ;;
	esac

	# search
	if test "$use_parallel" -ne 0
	then
		case "$arch" in
		alpha )
			AC_LANG_PUSH(C)
			CFLAGS="$test_cflags"
			CPPFLAGS="$test_cppflags"
			LDFLAGS="$c_ldflags"
			LIBS=""
			AC_SEARCH_LIBS(elan_baseInit, elan, mpi_libs="$LIBS")
			LIBS="-lmpi $mpi_libs"
			AC_SEARCH_LIBS(mpi_waitany_, fmpi, have_mpi=1
				       mpi_libs="$LIBS")
			AC_LANG_POP(C)

			LIBS="$mpi_libs"
			AC_SEARCH_LIBS(pmpi_init, pmpi, mpi_libs="$LIBS")
			;;
		esac

		if test "$have_mpi" -eq 0
		then
			# check for mpi
			unset ac_cv_search_mpi_init # clear cached value
			LDFLAGS="$test_ldflags"
			LIBS="$mpi_libs"
			AC_SEARCH_LIBS(mpi_init, mpi,
				       have_mpi=1 mpi_libs="$LIBS")
		fi
	fi
fi
echo setting MPI_LIBS... $mpi_libs

# final checks on available libraries
F77=$f90
FFLAGS="$test_fflags"
LDFLAGS="$test_ldflags"
parallel=0

if test "`echo $lapack_libs | grep mkl`" != ""
then
	mylib="lapack_mkl"
else
	LIBS="$lapack_libs $blas_libs"
	AC_SEARCH_LIBS(zggev, "", , mylib="lapack")
	AC_SEARCH_LIBS(dgemm, "", , mylib="blas_and_lapack")
fi

LIBS="$mass_libs"
AC_SEARCH_LIBS(vexp, "", try_dflags="$try_dflags -D__MASS")

if test "$arch" = "aix"
then
	LIBS="$fft_libs"
	AC_SEARCH_LIBS(dcft, "", have_fft=1)
fi

if test "$have_fft" -eq 0
then
	AC_LANG_PUSH(C)

	if test "$have_fft" -eq 0
	then
		try_dflags="$try_dflags -D__FFTW"

		CFLAGS="$test_cflags"
		LDFLAGS="$c_ldflags"
		LIBS="$fft_libs"
		AC_SEARCH_LIBS(fftwnd, "", have_fft=1 have_fftw=1)
	fi

	# find location of fftw.h (if needed)
	if test "$have_fftw" -ne 0
	then
		if test "$includefftw" = ""
		then
			# if libfftw is in /some/path/lib, then fftw.h is
			# probably in /some/path/include
			includefftw=`echo $fft_libs |
				sed 's/.*-L\([[^ ]]*\)  *-lfftw.*/\1/
				     s/lib/include/'`
		fi
		try_includedirs="$includefftw /usr/local/include \
				 /cineca/lib/fftw-2.1.3/fftw"
		for dir in none $try_includedirs
		do
			unset ac_cv_header_fftw_h # clear cached value
			if test "$dir" = "none"
			then
				try_ioption=
			else
				echo $ECHO_N "in $dir: " $ECHO_C
				try_ioption="-I$dir"
			fi
			CFLAGS="$test_cflags $try_ioption"
			CPPFLAGS="$test_cppflags $try_ioption"
			LDFLAGS="$c_ldflags"
			LIBS="$fft_libs"
			AC_CHECK_HEADER(fftw.h, have_fftw_h=1
					ifftw="$try_ioption")
			if test "$ac_cv_header_fftw_h" != "no"
			then break ; fi
		done
	fi

	# if no valid FFT library was found, use our implementation
	if test "$have_fft" -eq 0
	then
		try_dflags="$try_dflags -D__USE_INTERNAL_FFTW"
	fi

	if test "$have_fftw" -ne 0 && test "$have_fftw_h" -eq 0
	then
		try_dflags="$try_dflags -D__USE_INTERNAL_FFTW"
		AC_MSG_WARN([fftw library detected, but fftw.h not found])
	fi

	AC_LANG_POP(C)
fi

if test "$use_parallel" -ne 0
then
	F77=$mpif90
	LIBS="$mpi_libs"
	AC_SEARCH_LIBS(mpi_init, "", parallel=1
		       try_dflags="$try_dflags -D__MPI -D__PARA")
fi

AC_LANG_POP(Fortran 77)

if test "$dflags" = "" ; then dflags="$try_dflags" ; fi
echo setting DFLAGS... $dflags

# xlf compilers (AIX and powerpc) want comma-separated -D directives
if test "$xlf_flags" -ne 0
then
	fdflags="`echo $dflags | sed 's/  */,/g'`"
else
	fdflags="\$(DFLAGS)"
fi
echo setting FDFLAGS... $fdflags

if test "$ranlib" = "" ; then ranlib=echo ; fi
echo setting RANLIB... $ranlib
echo setting MYLIB... $mylib

# generate dependencies
if test -x $topdir/makedeps.sh
then
	echo checking dependencies...
	( cd $topdir ; ./makedeps.sh )
fi

# configure output messages
blas_line="BLAS_LIBS=$blas_libs"
lapack_line="LAPACK_LIBS=$lapack_libs"
fft_line="FFT_LIBS=$fft_libs"
if test "$have_fft" -ne 0 && test "$have_fftw" -ne 0 &&
	test "$have_fftw_h" -eq 0
then
	fft_line="$fft_line FFT_LIBS=$fft_libs\\
  WARNING: fftw library detected, but fftw.h not found"
fi
if test "$mpi_libs" != "" ; then mpi_line="MPI_LIBS=$mpi_libs"
else mpi_line="@delete@" ; fi
if test "$arch" = "aix" ; then mass_line="MASS_LIBS=$mass_libs"
else mass_line="@delete@" ; fi

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
AC_SUBST(mpicc)
AC_SUBST(cflags)
AC_SUBST(dflags)
AC_SUBST(fdflags)
AC_SUBST(cpp)
AC_SUBST(cppflags)
AC_SUBST(f90)
AC_SUBST(mpif90)
AC_SUBST(f90flags)
AC_SUBST(f90flags_noopt)
AC_SUBST(f77)
AC_SUBST(mpif77)
AC_SUBST(fflags)
AC_SUBST(fflags_noopt)
AC_SUBST(pre_fdflags)
AC_SUBST(imod)
AC_SUBST(ifftw)
AC_SUBST(ld)
AC_SUBST(ldflags)
AC_SUBST(blas_libs)
AC_SUBST(lapack_libs)
AC_SUBST(fft_libs)
AC_SUBST(mpi_libs)
AC_SUBST(mass_libs)
AC_SUBST(mylib)
AC_SUBST(ar)
AC_SUBST(arflags)
AC_SUBST(ranlib)
AC_SUBST(f90rule)
AC_SUBST(blas_line)
AC_SUBST(lapack_line)
AC_SUBST(fft_line)
AC_SUBST(mpi_line)
AC_SUBST(mass_line)
AC_SUBST(parallel_report)

AC_CONFIG_FILES(make.sys)
AC_CONFIG_FILES(configure.msg)

AC_OUTPUT

# final messages
sed '/@delete@/d' configure.msg
echo configure: success
