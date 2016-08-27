# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_CC], [

# candidate C and f77 compilers good for all cases
try_cc="cc gcc"

case "$arch:$f90_flavor" in
*:ifort* )
        try_cc="icc ecc $try_cc"
        ;;
*:pgf90 )
        try_cc="pgcc $try_cc"
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
        ;;
aix:*xlf* )
        try_cc="xlc $try_cc"
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
ppc64:*xlf* | ppc64-mn:*xlf* )
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

# compilation flags for all subsequent tests
test_cflags="`echo $cflags | sed 's/\$([[^)]]*)//g'`"

AC_SUBST(cflags)
])
