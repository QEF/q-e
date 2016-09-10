# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_BLAS], [

have_blas=0

# Flags for machine-specific libraries
have_acml=0
have_atlas=0
have_essl=0
have_mkl=0
  
AC_ARG_WITH(netlib,
   [AS_HELP_STRING([--with-netlib],
       [compile with Netlib LAPACK and BLAS (default: no)])],
    [if test "$withval" = "yes" ; then
      use_netlib=1
   else
      use_netlib=0
   fi],
   [use_netlib=0])

if test "$use_netlib" -eq 0
then
   if test "$blas_libs" = ""
   then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        case "$arch:$f90" in

        x86_64:path* | x86_64:openf95 | crayxt*:* )
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

        x86_64:pgf* )
                try_libdirs="/opt/acml*/pathscale64/lib/"
                try_libdirs="$ld_library_path $libdirs $try_libdirs"

                # Check first MKL...
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

                        # Check first MKL...
                        FFLAGS="$test_fflags"
                        LDFLAGS="$MKL_FLAGS $test_ldflags $try_loption"
                        LIBS=""

                        if test "$use_openmp" -eq 0; then
                              AC_SEARCH_LIBS(dgemm, mkl_intel_lp64,
                                 have_blas=1 have_mkl=1
                                 blas_libs="$try_loption $LIBS -lmkl_sequential -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core -ldl)
                        else
                              AC_SEARCH_LIBS(dgemm, mkl_intel_lp64,
                                 have_blas=1 have_mkl=1
                                 blas_libs="$try_loption $LIBS -lmkl_core -lmkl_pgi_thread"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core -ldl -lpthread -lm)
                        fi

                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi
                done

                # ... then ACML
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
                        LIBS=""
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
                        LIBS=""
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
                        LIBS=""
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
			LIBS=""
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

blas_line="BLAS_LIBS=$blas_libs" 
echo setting BLAS_LIBS... $blas_libs
  
AC_SUBST(blas_libs)
AC_SUBST(blas_line)
  
]
)
