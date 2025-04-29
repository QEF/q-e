# Copyright (C) 2001-2025 Quantum ESPRESSO Foundation

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

# pedantic flags implemented only for gcc
AC_ARG_ENABLE(pedantic,
   [AS_HELP_STRING([--enable-pedantic],
       [compile Fortran with pedantic flags (default: no)])],
   [if   test "$enableval" = "yes" ; then
      use_pedantic=1
   else
      use_pedantic=0
   fi],
   [use_pedantic=0])

# yet to be implemented
AC_ARG_ENABLE(shared,
   [AS_HELP_STRING([--enable-shared],
       [produce object files suitable for shared libraries (default: no)])],
   [if   test "$enableval" = "yes" ; then
      use_shared=1
   else
      use_shared=0
   fi],
   [use_shared=0])

# build static executables (implemented only for a few untested cases)
AC_ARG_ENABLE(static,
   [AS_HELP_STRING([--enable-static],
       [build static executables if possible (default: no)])],
   [if   test "$enableval" = "yes" ; then
      use_static=1
   else
      use_static=0
   fi],
   [use_static=0])

# check Fortran compiler flags
# have_cpp=0: use external C preprocessing for fortran code
# have_cpp=1: use C-like preprocessing in fortran compiler
have_cpp=1
xlf_flags=0

echo using F90... $f90

case "$arch:$f90_flavor" in
*:ifort* | *:ifx* )
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
        pre_fdflags="-fpp -allow nofpp_comments "
        ;;
arm:armflang )
        try_fflags="-O3 -mcpu=native $try_fflags"
        if test "$use_debug" -eq 1; then 
           try_fflags="$try_fflags -g" 
        fi   
        try_fflags_openmp="-fopenmp"
        try_f90flags="\$(FFLAGS) -cpp"
        try_ldflags="-g -mcpu=native"
        try_ldflags_openmp="-fopenmp"
        try_ldflags_static="-static -static-flang-libs"

        ;;
*:pgf* | *:nvfortran )
	try_fflags_nomain="-Mnomain"
        try_fflags="-fast"
        try_fflags_openmp="-mp"
        if test "$use_debug" -eq 1; then
           try_f90flags="-g -C -Ktrap=fp -Mcache_align -Mpreprocess -Mlarge_arrays"
        else
           try_f90flags="-fast -Mcache_align -Mpreprocess -Mlarge_arrays"
        fi
        try_foxflags="-fast -Mcache_align -Mpreprocess -Mlarge_arrays"
        try_fflags_noopt="-O0"
        try_ldflags=""
        try_ldflags_openmp="-mp"
        try_ldflags_static="-Bstatic"
        try_dflags="$try_dflags -D__PGI"
        have_cpp=1
        ;;
*:*gfortran )
	try_fflags="-O3 -g"
        if test "$use_debug" -eq 1; then
            try_fflags="-O3 -g  -Wall -fbounds-check -frange-check -finit-integer=987654321 -finit-real=nan -finit-logical=true -finit-character=64"
        fi
        if test "$use_pedantic" -eq 1; then
            try_fflags="-O2 -g -pedantic -Wall -Wextra -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all"
        fi
        if test "$f90_major_version" -ge "10"; then
 	   try_fflags="$try_fflags -fallow-argument-mismatch"
        fi
        try_fflags_openmp="-fopenmp"
        try_f90flags="\$(FFLAGS) -cpp"
        try_fflags_noopt="-O0 -g"
        try_ldflags="-g"
        try_ldflags_openmp="-pthread -fopenmp"
        try_ldflags_static="-static"
        ;;
*:flang )
        try_fflags="-O3"
        if test "$use_debug" -eq 1; then
            try_fflags="-O0 -g"
        fi
        try_fflags_nomain=""
        try_f90flags="\$(FFLAGS) -cpp"
        try_fflags_noopt="-O0 -g"
        try_dflags="$try_dflags -D_AOCC"
        try_ldflags=""
        ;;
# from now on: likely obsolete cases
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
craype*:cray* )
        try_fflags_nomain=""
        #NOTE: by default OpenMP is always ON (see crayftn man page)
        try_fflags_openmp="-homp"
        try_fflags="-O2"
        #NOTE: add '-rm' to get messages from crayftn about why
        #      optimizations have not been applied
        #      -x dir disable directives introduced by !DIR$
        try_f90flags="-eF -O3,fp3 -f free -x dir"
        try_fflags_noopt="-O0"
        try_ldflags_openmp="-homp"
        try_ldflags="-v"
        try_ldflags_static="-static"
        try_dflags="$try_dflags -D__CRAY"
        have_cpp=1
        ;;
craype*:pgf* )
# see comment above for pgf*
        try_fflags_nomain="-Mnomain"
        try_fflags_openmp="-mp"
        try_fflags="-O3"
        try_f90flags="-fast -Mcache_align -Mpreprocess -Mlarge_arrays"
        try_fflags_noopt="-O0"
        try_ldflags_openmp="-mp"
        try_ldflags="-v"
        try_dflags="$try_dflags -D__PGI"
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
        try_dflags="-D__XLF"
        pre_fdflags="-WF,"
        xlf_flags=1
        ;;
ppc64le:*xlf* )
    if test "$use_debug" -eq 1; then
        try_fflags="-g -C -qstrict -qdpc -qalias=nointptr -qarch=auto"
    else
        try_fflags="-O3 -qstrict -qdpc -qalias=nointptr -qarch=auto"
    fi
        try_fflags_openmp="-qsmp=noauto:omp"
        try_f90flags="\$(FFLAGS) -qsuffix=cpp=f90"
        try_fflags_noopt="-O0"
        try_ldflags=""
        try_ldflags_openmp="-qsmp=noauto:omp"
        try_dflags="-D__XLF"
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
        try_dflags="-D__XLF"
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
        try_dflags="-D__XLF"
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

if test "$use_static" -eq 1 ; then
  try_ldflags="$try_ldflags $try_ldflags_static" ; fi

# Flags are repeated, need better way to handle this ...
if test "$use_openmp" -eq 1 ; then
  try_f90flags="$try_f90flags $try_fflags_openmp"
  try_fflags="$try_fflags $try_fflags_openmp"
  try_ldflags="$try_ldflags $try_ldflags_openmp"
fi

if test "$fflags" = ""   ; then fflags=$try_fflags     ; fi
if test "$f90flags" = "" ; then f90flags=$try_f90flags ; fi
if test "try_foxflags" != ""; then foxflags=$try_foxflags; fi
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
AC_SUBST(foxflags)
])
