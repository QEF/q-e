AC_INIT(PWscf, 2.0, , pwscf)
AC_REVISION(0.01)
AC_PREREQ(2.53)

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

# configure for current directory by default
if test "$topdir" = "" ; then topdir="`pwd`" ; fi

# check system type (no cross-compilation for now)
AC_CANONICAL_BUILD

# check compiling environment
case $build in
*-pc-linux-* )
	is_linux=1
	try_f90="mpif90 ifort ifc pgf90 f90 mpif77 pgf77 f77 g77"
	try_ar="ar"
	try_arflags="ruv"
	try_dflags="-D__LINUX"
	;;
*-ibm-aix* )
	is_aix=1
	try_f90="mpxlf90_r mpxlf90 xlf90 f90 xlf f77 g77"
	try_ar="ar"
	try_arflags="-X64 ruv"
	try_dflags="-D__AIX"
	;;
mips-sgi-* )
	is_mips=1
	try_f90="f90 f77 g77"
	try_ar="ar"
	try_arflags="ruv"
	;;
sparc-sun-* )
	is_sparc=1
	try_f90="mpf90 f90 mpf77 f77 g77"
	try_ar="ar"
	try_arflags="ruv"
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

# check Fortran 90 compiler
if test "$f90" = "" ; then f90="$try_f90" ; fi
AC_PROG_F77($f90)
f90=$F77

case "$build:$f90" in
*-pc-linux-gnu:pgf* )
	is_pgi=1
	try_cc="pgcc gcc cc icc"
	try_f77="pgf77 $f90"
	;;
*-pc-linux-gnu:ifort | *-pc-linux-*:ifc )
	# extract compiler version
	echo $ECHO_N "checking version of $f90... $ECHO_C"
	ifc_version=`$f90 -V 2>&1 | grep "Intel(R)"`
	version=`echo $ifc_version | awk '{print $8}'`
	echo "${ECHO_T}$version"

	is_ifc=1
	try_cc="icc gcc cc mpicc pgcc"
	try_f77="$f90"
	;;
*-pc-linux-gnu:mpif* )
	# which compiler is actually under the hood?
	echo $ECHO_N "checking version of $f90... $ECHO_C"
	pgi_version=`$f90 -V 2>&1 | grep -e pgf`
	ifc_version=`$f90 -V 2>&1 | grep "Intel(R)"`
	if test "$pgi_version" != ""
	then
		name=`echo $pgi_version | awk '{print $1}'`
		version=`echo $pgi_version | awk '{print $2}'`
		echo "${ECHO_T}$name $version"

		is_pgi=1
		try_cc="mpicc pgcc gcc cc icc"
		try_f77="mpif77 pgf77 $f90"
	elif test "$ifc_version" != ""
	then
		version=`echo $ifc_version | awk '{print $8}'`
		echo "${ECHO_T}ifc $version"

		is_ifc=1
		try_cc="mpicc icc gcc cc pgcc"
		try_f77="mpif77 ifort ifc $f90"
	else
		echo "${ECHO_T}unknown"
		try_cc="mpicc gcc cc pgcc icc"
		try_f77="mpif77 $f90"
	fi
	;;
*-ibm-aix*:* )
	try_cc="mpcc_r mpcc xlc cc gcc"
	try_f77="mpxlf_r mpxlf xlf f77 $f90"
	;;
mips-sgi-*:* )
	try_cc="cc gcc"
	try_f77="f77 $f90"
	;;
sparc-sun-*:* )
	try_cc="mpcc cc gcc"
	try_f77="mpf77 f77 $f90"
	;;
* )
	# unknown
	try_cc="cc gcc"
	try_f77="f77 $f90"
	;;
esac

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

# check whether the C and Fortran compilers are compatible
case "$build:$cc:$f90" in
	# list supported combinations here
*-ibm-aix*:mpcc_r:mpxlf*_r | *-ibm-aix*:mpcc:mpxlf* | \
*-ibm-aix*:xlc:*xlf* | *-ibm-aix*:cc:*xlf* )            ;;
*-pc-linux-*:*cc:mpif90 )                               ;;
*-pc-linux-*:pgcc:pgf90 )                               ;;
*-pc-linux-*:icc:ifc | *-pc-linux-*:gcc:ifc | \
*-pc-linux-*:icc:ifort | *-pc-linux-*:gcc:ifort )       ;;
mips-sgi-*:cc:f90 )                                     ;;
# sparc-sun-*:mpcc:mpf90 | sparc-sun-*:cc:f90 )           ;;

* )
	AC_MSG_WARN(unsupported C/Fortran compilers combination)
	problems="$problems\\
		  unsupported C/Fortran compilers combination:\\
		  CC=$cc, F77=$f77, F90=$f90\\"
	;;
esac

# check compilers flags
case "$build:$f90" in
*-pc-linux-gnu:* )
	if test $is_pgi
	then
		try_fflags="-fast -r8"
		try_f90flags="\$(FFLAGS)"
		try_fflags_noopt="-O0"
		try_cpp="cpp"
		try_ldflags=""
		try_dflags="$try_dflags -D__PGI"
	elif test $is_ifc
	then
		# extract integer part of version
		if test `echo $version | sed 's/\..*//'` -lt 7
		then
			# old versions of ifc require this stuff
			echo $ECHO_N "setting up ifc environment... $ECHO_C"
			pcl_ph="D3"
			pcl_pw="PH PP Gamma PWNC PWCOND pwtools $pcl_ph"
			pcl_modules="PW FPMD CPV flib $pcl_pw"
			pcl_dot=". Modules upftools $pcl_modules"
			for dir in $pcl_dot
			do
				echo work.pc > $topdir/$dir/intel.pcl
			done
			for dir in $pcl_modules
			do
				echo ../Modules/work.pc \
				     >> $topdir/$dir/intel.pcl
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

			try_f90flags="\$(FFLAGS) -cl,./intel.pcl"
		else
			try_f90flags="\$(FFLAGS) -nomodule"
		fi
		try_fflags="-Vaxlib -O2 -tpp6"
		try_fflags_noopt="-O0"
		try_ldflags="-Vaxlib"
		try_dflags="$try_dflags -D__INTEL"
		pre_f_dflags="-fpp "
		have_cpp=1
	fi
	;;
*-ibm-aix*:mpxlf* | *-ibm-aix*:xlf* )
	try_fflags="-q64 -qalias=noaryovrlp -O3 -qstrict \
-qarch=auto -qtune=auto -qsuffix=cpp=f90 -qdpc -Q -qalias=nointptr"
	try_f90flags="\$(FFLAGS) -qfree=f90"
	try_fflags_noopt="-q64 -O0"
	try_ldflags="-q64"
	pre_f_dflags="-WF,"
	aix_flags=1
	have_cpp=1
	;;
mips-sgi-*:f90 )
	try_fflags="-mips4 -64 -O2 -r10000 -r8"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-mips4 -64 -O0"
	try_ldflags="-mips4 -64"
	try_dflags="-D__SGI -D__SGI64 -D__ORIGIN"
	pre_f_dflags="-cpp "
	have_cpp=1
	;;
sparc-sun-*:mpf90 | sparc-sun-*:f90 )
	try_fflags="-fast -dalign -xchip=ultra3 -xarch=v8plusb \
-xlic_lib=sunperf"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0"
	try_cpp="fpp"
	try_ldflags=""
	try_dflags="-DSUN"
	;;
*:g77 )
	try_fflags=""
	try_f90flags="\$(FFLAGS) -ff90 -ffree-form"
	try_fflags_noopt="-O0"
	try_ldflags=""
	;;
* )
	try_fflags="-O"
	try_f90flags="\$(FFLAGS)"
	try_fflags_noopt="-O0"
	try_ldflags=""
	;;
esac
if test "$fflags" = ""       ; then fflags=$try_fflags             ; fi
if test "$f90flags" = ""     ; then f90flags=$try_f90flags         ; fi
if test "$fflags_noopt" = "" ; then fflags_noopt=$try_fflags_noopt ; fi
echo setting FFLAGS... $fflags
echo setting F90FLAGS... $f90flags

case "$build:$cc" in
*-pc-linux-*:mpicc )
	if test $is_pgi
	then
		try_cflags="-fast"
		c_ldflags=""
	elif test $is_ifc
	then
		try_cflags="-O3"
		c_ldflags="-Vaxlib"
	fi
	;;
*-pc-linux-*:pgcc )
	try_cflags="-fast"
	c_ldflags=""
	;;
*-pc-linux-*:icc )
	try_cflags="-O3"
	c_ldflags="-Vaxlib"
	;;
*-ibm-aix*:mpcc_r | *-ibm-aix*:mpcc | *-ibm-aix*:xlc | *-ibm-aix*:cc )
	try_cflags="-q64 -O2"
	c_ldflags="-q64"
	pre_dflags="-Wp,"
	aix_flags=1
	;;
mips-sgi-*:cc )
	try_cflags="-mips4 -64 -O2 -r10000"
	;;
sparc-sun-*:mpcc | sparc-sun-*:cc )
	try_cflags="-fast -dalign -xchip=ultra3 -xarch=v8plusb \
-xlic_lib=sunperf"
	;;
*:gcc )
	try_cflags="-O3 -fomit-frame-pointer"
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
F77=$f90 # use Fortran 90 actually (autoconf knows Fortran 77 only!)

# search for libraries (unless LIBS has been set from the environment)
if test "$libs" = ""
then
	# build list of directories to search
	if test $is_aix
	then
		try_libdirs="/usr/local/lib /cineca/lib /cineca/lib/mass"
	elif test $is_linux && test $is_ifc
	then
		try_libdirs="/usr/local/lib
			     /opt/intel/mkl/mkl61/lib/32 /opt/intel/mkl/lib/32
			     /opt/intel/mkl61/lib/32
			     /cineca/lib /cineca/prod/intel/lib"
	elif test $is_linux && test $is_pgi
	then
		try_libdirs="/usr/local/lib
			     /usr/local/pgi/linux86/lib
			     /cineca/lib /cineca/prod/pgi/lib"
	elif test $is_mips
	then
		try_libdirs="/usr/local/lib /cineca/lib"
	fi

	# prepend directories specified by the user
	if test "$libdirs" != ""
	then
		try_libdirs="$libdirs $try_libdirs"
	fi

	if test "$LD_LIBRARY_PATH" != ""
	then
		ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`
		try_libdirs="$try_libdirs $ld_library_path"
	fi

	# use vendor libraries if available
	if test $is_aix
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
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS="$libs"
		AC_SEARCH_LIBS(dcft, essl, have_fft=1 have_essl=1 have_blas=1
			       libs="$LIBS")
	elif test $is_linux
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
$ac_cv_search_zggev -lmkl_ia32 -lguide $libs", , -lmkl_ia32 -lguide)
			if test "$ac_cv_search_zggev" != "no"
			then break ; fi
		done
	elif test $is_mips
	then
		# check for complib.sgimath
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS="$libs"
		AC_SEARCH_LIBS(dgemm, complib.sgimath, have_blas=1
			       libs="$LIBS")
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
				       libs="$try_loption $LIBS")
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
		# check for blas
		FFLAGS="$test_fflags"
		LDFLAGS="$test_ldflags"
		LIBS="$libs"
		AC_SEARCH_LIBS(dgemm, blas, have_blas=1 libs="$LIBS")
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

	# check for mpi
	LDFLAGS="$test_ldflags"
	LIBS="$libs"
	AC_SEARCH_LIBS(mpi_init, mpi, libs="$LIBS")
fi
echo setting LIBS... $libs

# final checks on available libraries
FFLAGS="$test_fflags"
LDFLAGS="$test_ldflags"
LIBS="$libs"

AC_SEARCH_LIBS(mpi_init, "", parallel=1
	       try_dflags="$try_dflags -D__MPI -D__PARA")

if test "`echo $libs | grep mkl`" != "" ; then mylib="lapack_mkl" ; fi
AC_SEARCH_LIBS(zggev, "", , mylib="blas_and_lapack")

AC_SEARCH_LIBS(zhegvx, "", try_dflags="$try_dflags -DHAS_ZHEGVX")
AC_SEARCH_LIBS(vexp, "", try_dflags="$try_dflags -D__MASS")

if test $is_aix
then
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

# AIX compilers want comma-separated -D directives
if test "$dflags" = "" ; then dflags="$try_dflags"                      ; fi
if test $aix_flags     ; then dflags="`echo $dflags | sed 's/  */,/g'`" ; fi
echo setting DFLAGS... $dflags
echo setting MYLIB... $mylib

# generate dependencies
if test -x $topdir/makedeps.sh
then
	echo $ECHO_N "checking module dependencies... $ECHO_C"
	( cd $topdir ; ./makedeps.sh )
	echo "${ECHO_T}done"
fi

# check sizes of types (not really needed)
CFLAGS="$test_cflags"
CPPFLAGS="$test_cppflags"
LDFLAGS="$c_ldflags"
LIBS=""
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(char)
AC_CHECK_SIZEOF(void *)

# export settings to generated files
AC_SUBST(topdir)
AC_SUBST(cc)
AC_SUBST(cflags)
AC_SUBST(dflags)
AC_SUBST(pre_dflags)
AC_SUBST(cpp)
AC_SUBST(cppflags)
AC_SUBST(f77)
AC_SUBST(f90)
AC_SUBST(fflags)
AC_SUBST(f90flags)
AC_SUBST(fflags_noopt)
AC_SUBST(pre_f_dflags)
AC_SUBST(ifftw)
AC_SUBST(ld)
AC_SUBST(ldflags)
AC_SUBST(libs)
AC_SUBST(mylib)
AC_SUBST(ar)
AC_SUBST(arflags)
AC_SUBST(f90rule)

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
echo For more information, read the INSTALL file or the PWscf manuals.
echo -----------------------------------------------------------------

if ! test $parallel
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
