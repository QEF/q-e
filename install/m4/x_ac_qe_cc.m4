# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_CC], [

# candidate C compilers good for all cases
try_cc="cc gcc"

case "$arch:$f90_flavor" in
*:ifort* )
        try_cc="icc $try_cc"
        ;;
*:pgf90 | *:nvfortran )
        try_cc="pgcc $try_cc"
        ;;
cray*:* )
        try_cc="cc"
        ;;
necsx:* )
        try_cc="sxcc"
        ;;
ppc64-bg*:*xlf90_r )
        try_cc="bgxlc_r"
        ;;
ppc64-bg*:*xlf90 )
        try_cc="bgxlc"
        ;;
ppc64:*xlf* | ppc64le:*xlf* )
        try_cc="xlc_r $try_cc"
        ;;
esac

# check serial C compiler
if test "$env_cc" = "" ; then cc="$try_cc" ; else cc="$env_cc"; fi
AC_PROG_CC($cc)
cc=$CC

echo setting CC... $cc

AC_SUBST(cc)

# tentative C and loader flags, good for many cases
try_cflags="-O3"
c_ldflags=""
try_cpp="cpp"

case "$arch:$cc" in
*:pgcc )
        try_cflags="-fast -Mpreprocess"
        # Workaround for BEEF compilation with PGI v.19 and previous 
        if test "$f90_flavor" = "pgf90"; then try_cflags="-c11 $try_cflags"; fi
        ;;
craype*:cc )
        # Actually we need something like is done for ftn to detect 
        # the proper compiler used (NdFilippo)
        try_cflags="-O3"
        ;;
necsx:* )
        #try_cflags="-D__SX6 \$(IFLAGS) \$(MODFLAGS)"
        try_cflags=""
        ;;
ppc64le:* )
        try_cflags="-O3"
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

# compilation flags for all subsequent tests
test_cflags="`echo $cflags | sed 's/\$([[^)]]*)//g'`"

AC_SUBST(cflags)
])
