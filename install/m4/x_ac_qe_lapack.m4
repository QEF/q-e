# Copyright (C) 2001-2025 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_LAPACK], [

if test "$have_mkl" -ne 0 || test "$have_armpl" -ne 0 || test "$have_aocl" -ne 0 || test "$have_essl" -ne 0
then
   # MKL or ARM libraries or AOCL or ESSL (obsolete?) found:
   # no need to check for lapack
   have_lapack=1
else
   # check for lapack
   have_lapack=0
fi
#
if test "$have_lapack" -eq 0
then
   if test "$lapack_libs" = ""
   then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        case "$arch:$f90" in

        necsx:* )
                # NECSX: OBSOLETE?
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
                # generic check for lapack (in several directories)
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
                done
        fi

   else
        # lapack provided in LAPACK_LIBS: not checked
	echo setting LAPACK from \$LAPACK_LIBS with no check ...  $lapack_libs
        have_lapack=1
   fi

fi

# No lapack library found: use internal lapack

if test "$have_lapack" -eq 0  ; then
    lapack_libs="\$(TOPDIR)/external/lapack/liblapack.a"
    echo setting LAPACK to internal library ...  $lapack_libs
    lapack_target="liblapack"
else
    lapack_target=""
fi
lapack_line="LAPACK_LIBS=$lapack_libs"

AC_SUBST(lapack_libs)
AC_SUBST(lapack_target)  
AC_SUBST(lapack_line)

AC_CONFIG_FILES(install/make_lapack.inc)
  
  ]
)
