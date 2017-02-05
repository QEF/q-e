# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_F90], [

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

# check Fortran compiler flags
# have_cpp=0: use external C preprocessing for fortran code
# have_cpp=1: use C-like preprocessing in fortran compiler
have_cpp=1
xlf_flags=0

echo using F90... $f90

case "$arch:$f90_flavor" in
ia32:ifort* | ia64:ifort* | x86_64:ifort* | mac686:ifort* | crayxt*:ifort* )
        try_fflags="-O2 -assume byterecl -g -traceback"
        if test "$use_debug" -eq 1; then
            try_fflags="$try_fflags -fpe0 -CB"
        fi
  	    try_fflags_nomain="-nofor_main"
        try_f90flags="\$(FFLAGS) -nomodule"
        try_fflags_noopt="-O0 -assume byterecl -g -traceback"
        try_ldflags=""
        try_ldflags_static="-static"
        if test "$f90_major_version" -ge "15"; then
            try_fflags_openmp="-qopenmp"
            try_ldflags_openmp="-qopenmp"
        else
            try_fflags_openmp="-openmp"
            try_ldflags_openmp="-openmp"
        fi
        pre_fdflags="-fpp "
        ;;
x86_64:nagfor* )
        try_fflags="-O3 -kind=byte -dcfuns -mismatch"
        if test "$use_debug" -eq 1; then
            try_fflags="$try_fflags -g"
        fi
        try_fflags_nomain=""
        try_fflags_openmp="-openmp"
        try_f90flags="-O3 -kind=byte -dcfuns -mismatch"
        try_fflags_noopt="-O0 -kind=byte -dcfuns -mismatch"
        try_ldflags=""
        try_ldflags_static="-unsharedrts"
        try_ldflags_openmp="-openmp"
        try_dflags="$try_dflags -D__NAG"
        have_cpp=0
        ;;
ia32:pgf* | ia64:pgf* | x86_64:pgf* )
	    try_fflags_nomain="-Mnomain"
        try_fflags="-fast -r8"
        try_fflags_openmp="-mp"
        try_f90flags="-fast -r8 -Mcache_align -Mpreprocess"
        try_fflags_noopt="-O0"
        try_ldflags=""
        try_ldflags_openmp="-mp"
        try_ldflags_static="-Bstatic"
        try_dflags="$try_dflags -D__PGI"
        have_cpp=1
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

# Flags are repeated, need better way to handle this ...
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

# compilation flags for all subsequent tests
# remove all $(...) because at least one compiler doesn't like them
# but if f90flags contains $(FFLAGS), substitute it
if test "`echo $f90flags | grep '$(FFLAGS)'`" != ""
then
        test_fflags="`echo $fflags $f90flags | sed 's/\$([[^)]]*)//g'`"
else
        test_fflags="`echo $f90flags | sed 's/\$([[^)]]*)//g'`"
fi

AC_SUBST(pre_fdflags)
AC_SUBST(f90flags)
AC_SUBST(fflags)
AC_SUBST(fflags_noopt)
AC_SUBST(fflags_nomain)
AC_SUBST(imod)

])
