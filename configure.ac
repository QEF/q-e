AC_INIT(PWscf, 2.1, , pwscf)
AC_REVISION(0.01)

# store variables from the environment (may be set or not)
topdir=$TOPDIR
cc=$CC
cpp=$CPP
cflags=$CFLAGS
cppflags=$CPPFLAGS
dflags=$DFLAGS
f77=$F77
f90=$F90
fflags=$FFLAGS
f90flags=$F90FLAGS
ld=$LD
ldflags=$LDFLAGS
libs=$LIBS
libdirs=$LIBDIRS
ar=$AR
arflags=$ARFLAGS
mylib=$MYLIB
includefftw=$INCLUDEFFTW

# command-line arguments
for arg in $*
do
	if   test "$arg" = "--disable-parallel" ; then no_parallel=1
	elif test "$arg" = "--disable-shared"   ; then no_shared=1
	elif test "$arg" = "--enable-shared"    ; then yes_shared=1
	fi
done

# configure for current directory by default
if test "$topdir" = "" ; then topdir="`pwd`" ; fi

# check system type (no cross-compilation for now)
AC_CANONICAL_BUILD

# check compiling environment
case $build in
ia64-*-linux-gnu )
	is_linux64=1
	try_f90_parallel="mpif90"
	try_f90="ifort ifc efc f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__LINUX64"
	;;
x86_64-*-linux-gnu )
	is_linux64=1
	try_f90_parallel="mpif90"
	try_f90="pgf90 f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__LINUX64"
	;;
*-pc-linux-gnu )
	is_linux=1
	try_f90_parallel="mpif90"
	try_f90="ifort ifc pgf90 g95 f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__LINUX"
	# use static libraries unless explicitly asked
	if test ! $yes_shared ; then no_shared=1 ; fi
	;;
*-ibm-aix* )
	is_aix=1
	try_f90_parallel="mpxlf90_r mpxlf90"
	try_f90="xlf90_r xlf90 f90"
	try_ar="ar"
	try_arflags="-X64 ruv"
	try_dflags="-D__AIX"
	;;
mips-sgi-* )
	is_mips=1
	try_f90="f90"
	try_ar="ar"
	try_arflags="ruv"
	;;
alphaev*-dec-osf* )
	is_alpha=1
	try_f90="f90"
	try_ar="ar"
	try_arflags="ruv"
	;;
sparc-sun-* )
	is_sparc=1
	try_f90_parallel="mpf90"
	try_f90="f90"
	try_ar="ar"
	try_arflags="ruv"
	;;
*cray-unicosmp* )
	is_x1=1
	try_f90="ftn"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__X1"
	;;
powerpc-apple-darwin* )
	is_mac=1
	try_f90_parallel="mpxlf90 mpf90"
	try_f90="xlf90 f90"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__MAC"
	ranlib="ranlib"
	;;
* )
	AC_MSG_WARN(unsupported architecture)
	problems="$problems\\
		  unsupported architecture:\\
		  $build\\"
	try_f90="f90"
	try_ar="ar"
	try_arflags="ruv"
	;;
esac
if ! test $no_parallel ; then try_f90="$try_f90_parallel $try_f90" ; fi

# check Fortran 90 compiler
if test "$f90" = "" ; then f90="$try_f90" ; fi
unset F77 # otherwise it may interfere with the check
AC_PROG_F77($f90)
f90=$F77
f90_id=$f90 # this should be the compiler's brand name

case "$build:$f90" in
ia64-*-linux-gnu:* | x86_64-*-linux-gnu:* | *-pc-linux-gnu:* )
	# double-check compiler
	# because on GNU/Linux systems they often have nonstandard names
	echo $ECHO_N "checking version of $f90... $ECHO_C"
	ifort_version=`$f90 -V 2>&1 | grep "Intel(R)"`
	pgf_version=`$f90 -V 2>&1 | grep "^pgf"`
	g95_version=`$f90 -v 2>&1 | grep "version"`
	if test "$ifort_version" != ""
	then
		version=`echo $ifort_version | awk '{print $8}'`
		if test `echo $version | sed 's/\..*//'` -lt 8
		then
			f90_id=ifc
		else
			f90_id=ifort
		fi
		echo "${ECHO_T}$f90_id $version"

		is_ifort=1
		try_cc="icc gcc cc"
		try_f77="ifort ifc efc $f90"
	elif test "$pgf_version" != ""
	then
		f90_id=`echo $pgf_version | awk '{print $1}'`
		version=`echo $pgf_version | awk '{print $2}'`
		echo "${ECHO_T}$f90_id $version"

		is_pgf=1
		try_cc="pgcc gcc cc"
		try_f77="pgf77 $f90"
	elif test "$g95_version" != ""
	then
		f90_id=g95
		version=`echo $g95_version | awk '{print $3}'`
		echo "${ECHO_T}g95 $version"

		is_g95=1
		try_cc="gcc cc"
		try_f77="$f90"
	else
		echo "${ECHO_T}unknown"
		try_cc="gcc cc"
		try_f77="f77 $f90"
	fi
	;;
*-ibm-aix*:*xlf*_r )
	try_cc="xlc_r cc gcc"
	try_f77="xlf_r f77 $f90"
	;;
*-ibm-aix*:*xlf* )
	try_cc="xlc cc gcc"
	try_f77="xlf f77 $f90"
	;;
mips-sgi-*:* )
	try_cc="cc gcc"
	try_f77="f77 $f90"
	;;
alphaev*-dec-osf* )
	try_cc="cc gcc"
	try_f77="f77 $f90"
	;;
sparc-sun-*:* )
	try_cc="cc gcc"
	try_f77="f77 $f90"
	;;
*cray-unicosmp*:ftn )
	try_cc="cc"
	try_f77="$f90"
	;;
powerpc-apple-darwin*:xlf90 )
	try_cc="gcc"
	try_f77="xlf f77 $f90"
	;;
powerpc-apple-darwin* )
	try_cc="gcc"
	try_f77="f77 $f90"
	;;
* )
	# unknown, try these
	try_cc="cc gcc"
	try_f77="f77 $f90"
	;;
esac

case "$build:$f90" in
*:mpif90 )
	try_cc_parallel="mpicc"
	try_f77_parallel="mpif77"
	;;
*-ibm-aix*:mpxlf*_r )
	try_cc_parallel="mpcc_r"
	try_f77_parallel="mpxlf_r"
	;;
*-ibm-aix*:mpxlf* )
	try_cc_parallel="mpcc"
	try_f77_parallel="mpxlf"
	;;
sparc-sun-*:mpf90 )
	try_cc_parallel="mpcc"
	try_f77_parallel="mpf77"
	;;
powerpc-apple-darwin*:mp* )
	try_cc_parallel="mpxlc mpcc"
	try_f77_parallel="mpxlf77 mpf77"
	;;
esac
if ! test $no_parallel
then
	try_f77="$try_f77_parallel $try_f77"
	try_cc="$try_cc_parallel $try_cc"
fi

# clear cached values
unset F77 ac_cv_prog_ac_ct_F77 ac_cv_f77_compiler_gnu ac_cv_prog_f77_g

# Fortran 77 compiler
if test "$f77" = "" ; then f77="$try_f77" ; fi
AC_PROG_F77($f77)
f77=$F77

# C compiler
if test "$cc" = "" ; then cc="$try_cc" ; fi
AC_PROG_CC($cc)
cc=$CC
cc_id=$cc # this should be the compiler's brand name

case "$build:$cc" in
ia64-*-linux-gnu:* | x86_64-*-linux-gnu | *-pc-linux-gnu:* )
	# double-check compiler
	# because on GNU/Linux systems they often have nonstandard names
	echo $ECHO_N "checking version of $cc... $ECHO_C"
	icc_version=`$cc -V 2>&1 | grep "Intel(R)"`
	pgcc_version=`$cc -V 2>&1 | grep "^pgcc"`
	gcc_version=`$cc -v 2>&1 | grep "version"`
	if test "$icc_version" != ""
	then
		cc_id=icc
		cversion=`echo $icc_version | awk '{print $8}'`
		echo "${ECHO_T}icc $cversion"
	elif test "$pgcc_version" != ""
	then
		cc_id=`echo $pgcc_version | awk '{print $1}'`
		cversion=`echo $pgcc_version | awk '{print $2}'`
		echo "${ECHO_T}$cc_id $cversion"
	elif test "$gcc_version" != ""
	then
		cc_id=gcc
		cversion=`echo $gcc_version | awk '{print $3}'`
		echo "${ECHO_T}gcc $cversion"
	else
		echo "${ECHO_T}unknown"
	fi
	;;
esac

# check whether the C and Fortran compilers are compatible
case "$build:$cc_id:$f90_id" in
	# list supported combinations here
ia64-*-linux-gnu:icc:ifort )                            ;;
ia64-*-linux-gnu:gcc:ifort | ia64-*-linux-gnu:gcc:efc ) ;;
x86_64-*-linux-gnu:pgcc:pgf90 )                         ;;
*-pc-linux-gnu:icc:ifort | *-pc-linux-*:icc:ifc )       ;;
*-pc-linux-gnu:gcc:ifort | *-pc-linux-*:gcc:ifc )       ;;
*-pc-linux-gnu:pgcc:pgf90 )                             ;;
*-pc-linux-gnu:gcc:g95 )                                ;;
*-ibm-aix*:mpcc*:mpxlf* )                               ;;
*-ibm-aix*:xlc*:*xlf* | *-ibm-aix*:cc:*xlf* )           ;;
mips-sgi-*:cc:f90 )                                     ;;
alphaev*-dec-osf*:cc:f90 )                              ;;
powerpc-apple-darwin*:mpcc:mpf90 )                      ;;
powerpc-apple-darwin*:gcc:xlf90 )                       ;;
*cray-unicosmp*:cc:ftn )                                ;;
# sparc-sun-*:mpcc:mpf90 | sparc-sun-*:cc:f90 )           ;;

* )
	AC_MSG_WARN(unsupported C/Fortran compilers combination)
	problems="$problems\\
		  unsupported C/Fortran compilers combination:\\
		  CC=$cc, F77=$f77, F90=$f90\\"
	;;
esac

# check Fortran compiler flags
case "$build:$f90_id" in
ia64-*-linux-gnu:ifort | ia64-*-linux-gnu:efc )
	try_fflags="-Vaxlib -O2 -assume byterecl"
	try_f90flags="\$(FFLAGS) -nomodule"
	try_fflags_noopt="-O0"
	try_ldflags="-Vaxlib"
	try_ldflags_static="-static"
	try_dflags="$try_dflags -D__INTEL"
	pre_fdflags="-fpp "
	have_cpp=1
	;;
x86_64-*-linux-gnu:pgf* )
	try_fflags="-fast -r8"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0"
	try_ldflags=""
	try_ldflags_static="-Bstatic"
	try_dflags="$try_dflags -D__PGI"
	;;
*-pc-linux-gnu:ifort )
        try_fflags="-Vaxlib -O2 -tpp6 -assume byterecl"
	try_f90flags="\$(FFLAGS) -nomodule"
	try_fflags_noopt="-O0"
	try_ldflags="-Vaxlib"
	try_ldflags_static="-static"
	try_dflags="$try_dflags -D__INTEL"
	pre_fdflags="-fpp "
	have_cpp=1
	;;
*-pc-linux-gnu:ifc )
        try_fflags="-Vaxlib -O2 -tpp6"
	try_f90flags="\$(FFLAGS) -nomodule"
	try_fflags_noopt="-O0"
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
		pcl_ph="D3"
		pcl_pw="PH PP Gamma PWNC PWCOND pwtools $pcl_ph"
		pcl_modules="PW FPMD CPV flib upftools $pcl_pw"
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
*-pc-linux-gnu:pgf* )
	try_fflags="-fast -r8"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0"
	try_ldflags=""
	try_ldflags_static="-Bstatic"
	try_dflags="$try_dflags -D__PGI"
	;;
*-pc-linux-gnu:g95 )
	try_fflags="-O3"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0"
	try_ldflags=""
	#try_ldflags_static="-static"
	try_dflags="$try_dflags -D__G95"
	;;
*-ibm-aix*:*xlf* )
	try_fflags="-q64 -qalias=noaryovrlp -O3 -qstrict \
-qarch=auto -qtune=auto -qsuffix=cpp=f90 -qdpc -Q -qalias=nointptr"
	try_f90flags="\$(FFLAGS) -qfree=f90"
	try_fflags_noopt="-q64 -O0"
	try_ldflags="-q64"
	# try_ldflags_static="-bstatic"
	pre_fdflags="-WF,"
	xlf_flags=1
	have_cpp=1
	;;
mips-sgi-*:f90 )
	try_fflags="-mips4 -64 -O2 -r10000 -r8"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-mips4 -64 -O0"
	try_ldflags="-mips4 -64"
	try_dflags="-D__SGI -D__SGI64 -D__ORIGIN"
	pre_fdflags="-cpp "
	have_cpp=1
	;;
alphaev*-dec-osf*:f90 )
	have_cpp=1
	try_dflags="-D__ALPHA"
	pre_fdflags="-cpp "
	try_fflags="-O -real_size 64 -align dcommons -align records"
	try_f77flags="\$(FFLAGS)"
	try_f90flags="\$(FFLAGS) -free"
	;;
sparc-sun-*:mpf90 | sparc-sun-*:f90 )
	try_fflags="-fast -dalign -xchip=ultra3 -xarch=v8plusb \
-xlic_lib=sunperf"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0"
	try_ldflags=""
	try_dflags="-DSUN"
	;;
*cray-unicosmp*:ftn )
        try_fflags="-s default64 -dp -rma -e0 -O scalar3,stream3,vector3,nointerchange -Z -O inline5"
        try_f90flags="\$(FFLAGS)"
        try_fflags_noopt="FFLAGS = -s default64 -dp -rma -e0 -O0 -Z"
        try_ldflags="-s default64 -f nan64"
        try_dflags="$try_dflags -D__X1_COA"
        pre_fdflags="-e Z -F"
        have_cpp=1
	;;
powerpc-apple-darwin*:* )
	try_fflags="-O4 -qarch=auto -qtune=auto -qsuffix=cpp=f90 -qdpc -qalias=nointptr"
	try_f90flags="\$(FFLAGS) -qfree=f90"
	try_fflags_noopt="-O0"
	try_ldflags=""
	pre_fdflags="-WF,"
	xlf_flags=1
	have_cpp=1
	;;	
*:g77 )
	# g77 is known not to work, but one might want to try anyway
	try_fflags=""
	try_f90flags="\$(FFLAGS) -ff90 -ffree-form"
	try_fflags_noopt="-O0"
	try_ldflags=""
	;;
* )
	# unknown, try these
	try_fflags="-O"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0"
	try_ldflags=""
	;;
esac
if test $no_shared ; then try_ldflags="$try_ldflags $try_ldflags_static" ; fi

if test "$fflags" = ""       ; then fflags=$try_fflags             ; fi
if test "$f90flags" = ""     ; then f90flags=$try_f90flags         ; fi
if test "$fflags_noopt" = "" ; then fflags_noopt=$try_fflags_noopt ; fi
echo setting FFLAGS... $fflags
echo setting F90FLAGS... $f90flags

case "$build:$cc_id" in
ia64-*-linux-gnu:icc )
	try_cflags="-O3"
	c_ldflags="-Vaxlib"
	;;
x86_64-*-linux-gnu:pgcc )
	try_cflags="-fast"
	c_ldflags=""
	try_cpp="cpp"
	;;
*-pc-linux-gnu:icc )
	try_cflags="-O3"
	c_ldflags="-Vaxlib"
	;;
*-pc-linux-gnu:pgcc )
	try_cflags="-fast"
	c_ldflags=""
	try_cpp="cpp"
	;;
*-pc-linux-gnu:gcc )
	try_cflags="-O3 -fomit-frame-pointer"
	c_ldflags=""
	try_cpp="cpp"
	;;
*-ibm-aix*:mpcc* | *-ibm-aix*:xlc* | *-ibm-aix*:cc )
	try_cflags="-q64 -O2"
	c_ldflags="-q64"
	;;
mips-sgi-*:cc )
	try_cflags="-mips4 -64 -O2 -r10000"
	;;
alphaev*-dec-osf*:cc )
	try_cflags="-O"
	;;
sparc-sun-*:mpcc | sparc-sun-*:cc )
	try_cflags="-fast -dalign -xchip=ultra3 -xarch=v8plusb \
-xlic_lib=sunperf"
	try_cpp="fpp"
	;;
*cray-unicosmp*:cc )
        try_cflags=""
        c_ldflags=""
        ;;
powerpc-apple-darwin*:gcc )
	try_cflags="-O3 -fomit-frame-pointer -I/usr/include/malloc"
	try_cpp="cpp"
	;;
powerpc-apple-darwin*:* )
	try_cflags="-O4"
	;;
*:gcc )
	try_cflags="-O3 -fomit-frame-pointer"
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

# linker and archiver (no tests)
if test "$ld" = ""       ; then ld="$f90"              ; fi
if test "$ldflags" = ""  ; then ldflags="$try_ldflags" ; fi
if test "$ar" = ""       ; then ar="$try_ar"           ; fi
if test "$arflags" = ""  ; then arflags="$try_arflags" ; fi
echo setting LD... $ld
echo setting LDFLAGS... $ldflags
echo setting AR... $ar
echo setting ARFLAGS... $arflags

# compilation rules
echo $ECHO_N "checking whether Fortran files must be preprocessed... $ECHO_C"
if test $have_cpp
then
	f90rule="\$(F90) \$(F90FLAGS) \$(MODULEFLAG) -c \$<"
	echo "${ECHO_T}no"
else
	f90rule="\$(CPP) \$(CPPFLAGS) \$*.f90 \$*.F90 \\
	\$(F90) \$(F90FLAGS) \$(MODULEFLAG) -c \$*.F90 -o \$*.o"
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

# search for libraries (unless LIBS has been set from the environment)
if test "$libs" = ""
then
	# build list of directories to search
	if test $is_linux64 && test $is_ifort
	then
		try_libdirs="/usr/local/lib
			     /opt/intel/mkl70/lib/64
			     /opt/intel/mkl/mkl61/lib/64 /opt/intel/mkl/lib/64
			     /opt/intel/mkl61/lib/64"
	elif test $is_linux64 && test $is_pgf
	then
		try_libdirs="/usr/local/lib
			     /usr/local/pgi/linux86/lib
			     /opt/pgi/linux86/5.1/lib
			     /cineca/lib /cineca/prod/pgi/lib"
	elif test $is_linux && test $is_ifort
	then
		try_libdirs="/usr/local/lib
			     /opt/intel/mkl70/lib/32
			     /opt/intel/mkl/mkl61/lib/32 /opt/intel/mkl/lib/32
			     /opt/intel/mkl61/lib/32
			     /cineca/lib /cineca/prod/intel/lib"
	elif test $is_linux && test $is_pgf
	then
		try_libdirs="/usr/local/lib
			     /usr/local/pgi/linux86/lib
			     /opt/pgi/linux86/5.1/lib
			     /cineca/lib /cineca/prod/pgi/lib"
	elif test $is_aix
	then
		try_libdirs="/usr/local/lib /cineca/lib /cineca/lib/mass"
	elif test $is_mips
	then
		try_libdirs="/usr/local/lib /cineca/lib"
	elif test $is_mac
	then
		try_libdirs="/usr/local/lib /sw/lib"
	else
		try_libdirs="/usr/local/lib"
	fi

	# prepend directories specified by the user
	if test "$libdirs" != ""
	then
		try_libdirs="$libdirs $try_libdirs"
	fi

	# check directories in LD_LIBRARY_PATH too
	# (maybe they are already searched by default, but I'm not sure)
	if test "$LD_LIBRARY_PATH" != ""
	then
		ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`
		try_libdirs="$try_libdirs $ld_library_path"
	fi

	# use vendor libraries if available
	if test $is_linux64 && test $is_ifort
	then
		# check for mkl (in several directories)
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
			LIBS="$libs"
			AC_SEARCH_LIBS(zggev, mkl_lapack,
				       have_blas=1 have_lapack=1
				       libs="$try_loption \
$ac_cv_search_zggev -lmkl_ipf -lguide $libs", , -lmkl_ipf -lguide)
			if test "$ac_cv_search_zggev" != "no"
			then break ; fi
		done
	elif test $is_linux && test $is_ifort
	then
		# check for mkl (in several directories)
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
			LIBS="$libs"
			AC_SEARCH_LIBS(zggev, mkl_lapack,
				       have_blas=1 have_lapack=1
				       libs="$try_loption \
$ac_cv_search_zggev -lmkl_ia32 -lguide -lpthread $libs", ,
				       -lmkl_ia32 -lguide -lpthread)
			if test "$ac_cv_search_zggev" != "no"
			then break ; fi
		done
	elif test $is_aix
	then
		# check for mass (in several directories)
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
			LIBS="$libs"
			AC_SEARCH_LIBS(vexp, massvp4 massv,
				       libs="$try_loption $ac_cv_search_vexp \
-lmass $libs", , -lmass)
			if test "$ac_cv_search_vexp" != "no" ; then break ; fi
		done

		# check for essl
		unset ac_cv_search_dcft # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS="$libs"
		AC_SEARCH_LIBS(dcft, essl, have_fft=1 have_essl=1 have_blas=1
			       libs="$LIBS")
	elif test $is_mips
	then
		# check for complib.sgimath
		unset ac_cv_search_dgemm # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS="$libs"
		AC_SEARCH_LIBS(dgemm, complib.sgimath, have_blas=1
			       libs="$LIBS")
	elif test $is_alpha
	then
		# check for cxml (old dxml)
		unset ac_cv_search_dgemm # clear cached value
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS="$libs"
		AC_SEARCH_LIBS(dgemm, cxml, have_blas=1 have_lapack=1 
			       libs="$LIBS")

		# check for parallel libraries
		if ! test $no_parallel
		then
			AC_LANG_PUSH(C)
			CFLAGS="$test_cflags"
			CPPFLAGS="$test_cppflags"
			LDFLAGS="$c_ldflags"
			LIBS="$libs"
			AC_SEARCH_LIBS(elan_baseInit, elan, libs="$LIBS")
			AC_SEARCH_LIBS(mpi_waitany_, fmpi, have_mpi=1
				       libs="$ac_cv_search_mpi_waitany_ \
-lmpi $libs", , -lmpi)
			AC_LANG_POP(C)

			LIBS="$libs"
			AC_SEARCH_LIBS(pmpi_init, pmpi, libs="$LIBS")
		fi
	fi

	if ! test $have_fft
	then
		AC_LANG_PUSH(C)

		# check for fftw (in several directories)
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
			LIBS="$libs"
			AC_SEARCH_LIBS(fftwnd, fftw, have_fftw=1
				       libs="$try_loption $LIBS", , -lm)
			if test "$ac_cv_search_fftwnd" != "no"
			then break ; fi
		done

		AC_LANG_POP(C)
	fi

	if ! test $have_blas
	then
		# check for atlas (in several directories)
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
			LIBS="$libs"
			AC_SEARCH_LIBS(dgemm, f77blas, have_blas=1
				       libs="$try_loption \
$ac_cv_search_dgemm -latlas $libs", , -latlas)
			if test "$ac_cv_search_dgemm" != "no"
			then break ; fi
		done
	fi

	if ! test $have_blas
	then
		# check for blas (in several directories)
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
			LIBS="$libs"
			AC_SEARCH_LIBS(dgemm, blas, have_blas=1
				       libs="$try_loption $LIBS")
			if test "$ac_cv_search_dgemm" != "no"
			then break ; fi
		done
	fi

	if ! test $have_lapack
	then
		# check for lapack (in several directories)
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
			LIBS="$libs"
			AC_SEARCH_LIBS(zggev, lapack, have_lapack=1
				       libs="$try_loption $LIBS")
			if test "$ac_cv_search_zggev" != "no"
			then break ; fi
		done
	fi
	# essl must precede lapack (if present)
	if test $have_essl && test $have_lapack
	then libs="-lessl $libs" ; fi

	if ! test $no_parallel && ! test $have_mpi
	then
		# check for mpi
		unset ac_cv_search_mpi_init # clear cached value
		LDFLAGS="$test_ldflags"
		LIBS="$libs"
		AC_SEARCH_LIBS(mpi_init, mpi, have_mpi=1 libs="$LIBS")
	fi
fi
echo setting LIBS... $libs

# final checks on available libraries
FFLAGS="$test_fflags"
LDFLAGS="$test_ldflags"
LIBS="$libs"

if ! test $no_parallel
then
	AC_SEARCH_LIBS(mpi_init, "", parallel=1
		       try_dflags="$try_dflags -D__MPI -D__PARA")
fi

if test "`echo $libs | grep mkl`" != "" ; then mylib="lapack_mkl" ; fi
AC_SEARCH_LIBS(zggev, "", , mylib="blas_and_lapack")

if test $is_aix
then
	AC_SEARCH_LIBS(vexp, "", try_dflags="$try_dflags -D__MASS")
	AC_SEARCH_LIBS(dcft, "", have_fft=1)
fi

AC_LANG_PUSH(C)

if ! test $have_fft
then
	try_dflags="$try_dflags -D__FFTW"

	CFLAGS="$test_cflags"
	LDFLAGS="$c_ldflags"
	AC_SEARCH_LIBS(fftwnd, "", have_fft=1 have_fftw=1)
fi

# find location of fftw.h (if needed)
if test $have_fftw
then
	if test "$includefftw" = ""
	then
		# if libfftw is in /some/path/lib, then fftw.h is
		# likely to be in /some/path/include
		includefftw=`echo $libs |
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
		LIBS="$libs"
		AC_CHECK_HEADER(fftw.h, have_fftw_h=1 ifftw="$try_ioption")
		if test "$ac_cv_header_fftw_h" != "no" ; then break ; fi
	done
fi

# if no valid FFT library was found, use our implementation
if ! test $have_fft
then
	try_dflags="$try_dflags -D__USE_INTERNAL_FFTW"
fi

if test $have_fftw && ! test $have_fftw_h
then
	try_dflags="$try_dflags -D__USE_INTERNAL_FFTW"
	AC_MSG_WARN([fftw library detected, but fftw.h not found])
fi

AC_LANG_POP(C)
AC_LANG_POP(Fortran 77)

if test "$dflags" = "" ; then dflags="$try_dflags" ; fi
echo setting DFLAGS... $dflags

# xlf compilers (AIX and powerpc) want comma-separated -D directives
if test $xlf_flags
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
	echo $ECHO_N "checking module dependencies... $ECHO_C"
	( cd $topdir ; ./makedeps.sh )
	echo "${ECHO_T}done"
fi

# export settings to generated files
AC_SUBST(cc)
AC_SUBST(cflags)
AC_SUBST(dflags)
AC_SUBST(fdflags)
AC_SUBST(cpp)
AC_SUBST(cppflags)
AC_SUBST(f77)
AC_SUBST(f90)
AC_SUBST(fflags)
AC_SUBST(f90flags)
AC_SUBST(fflags_noopt)
AC_SUBST(pre_fdflags)
AC_SUBST(ifftw)
AC_SUBST(ld)
AC_SUBST(ldflags)
AC_SUBST(libs)
AC_SUBST(mylib)
AC_SUBST(ar)
AC_SUBST(arflags)
AC_SUBST(f90rule)
AC_SUBST(ranlib)

AC_CONFIG_FILES(make.sys make.rules)

AC_OUTPUT

# final warnings
echo -----------------------------------------------------------------
echo PWscf can take advantage of several optimized numerical libraries
echo \(essl, fftw, mkl...\).  This configure script attempts to find them,
echo but may fail if they have been installed in non-standard locations.
echo
echo The following libraries have been found:
echo "  LIBS=$libs"

if test $have_fft && test $have_fftw && ! test $have_fftw_h
then
	echo "  WARNING: fftw library detected, but fftw.h not found"
fi

echo
echo If any libraries are missing, you may specify a list of directories
echo to search and retry, as follows:
echo "  ./configure LIBDIRS=\"list of directories, separated by spaces\""
echo
echo For more information, read the README.INSTALL file or the PWscf manuals.
echo -----------------------------------------------------------------

if ! test $no_parallel && ! test $parallel
then
	echo ----------------------------------------------
	echo WARNING: parallel environment not detected
	echo this program will run in single-processor mode
	echo ----------------------------------------------
fi

if test "$problems" != ""
then
	echo ----------------------------------------------
	echo WARNING: the following problems were detected:
	echo $problems | tr \\\\ \\n
	echo you may not be able to compile this program
	echo ----------------------------------------------
fi
