# Copyright (C) 2001-2014 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_LAPACK], [

have_lapack=0

AC_ARG_WITH(internal-lapack,
   [AS_HELP_STRING([--with-internal-lapack],
       [compile with internal lapack (default: no)])],
    [if   test "$withval" = "yes" ; then
      use_internal_lapack=1
   else
      use_internal_lapack=0
   fi],
   [use_internal_lapack=0])
   
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

  lapack_line="LAPACK_LIBS=$lapack_libs"

  AC_SUBST(lapack_libs)
  AC_SUBST(lapack_libs_switch)  
  AC_SUBST(lapack_line)
  
  AC_CONFIG_FILES(install/make_lapack.inc)
  
  ]
)
