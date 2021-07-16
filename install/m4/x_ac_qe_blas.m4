# Copyright (C) 2001-2020 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_BLAS], [

have_blas=0

# Flags for machine-specific libraries
have_acml=0
have_atlas=0
have_essl=0
have_mkl=0
have_armpl=0 

if test "$blas_libs" != ""
then
    echo setting BLAS from \$BLAS_LIBS with no check ...  $blas_libs
    have_blas=1
else
    # check directories in LD_LIBRARY_PATH too
    # (maybe they are already searched by default: useless?)
    ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

    case "$arch:$f90" in

    # search for architecture-specific libraries
    
    x86_64:* | mac686:* )
            #
            # search for MKL in directory $MKL_ROOT
	    #
	    # Following architectures no longer supported:
	    #   ia64  $MKLROOT/lib/64   -lmkl_gf_ipf, -lmkl_intel_ipf
	    #   ia32  $MKLROOT/lib/ia32 -lmkl_gf    , -lmkl_intel
            #
            if test "$MKLROOT" == ""; then
               MKLROOT=/opt/intel/mkl
            fi
	    case "$f90" in
	       ifort* )
      		    mkl_lib="mkl_intel_lp64"
      		    mkl_omp="mkl_intel_thread"
		    if test "$arch" == "mac686"; then
		       add_mkl_flag="-openmp"
		       add_mkl_lib="-lpthread"
		       add_mkl_omp="-lpthread"
		    fi
		    ;;
	       gfortran* )
      		    mkl_lib="mkl_gf_lp64"
      		    mkl_omp="mkl_gnu_thread"
		    ;;
	       nvfortran* )
                    # NB: next two can be replaced by flag "-Mmkl"
      		    mkl_lib="mkl_intel_lp64"
                    mkl_omp="mkl_intel_thread"
		    # NB: with nvidia hpc sdk 2020, linking to threaded mkl
		    # v.19.1 update 4 fails due to a missing symbol
		    ;;
	       pgf* )
                    # For obsolete PGI versions (superseded by nvfortran)
                    pgf_version=`$mpif90 -V 2>&1 | sed '/^$/d' | grep "^pgf" | cut -d ' ' -f2`
                    # From version 19.1, the new llvm backend requires linking to mkl_intel_thread
                    ompimp=""
                    AS_VERSION_COMPARE([$pgf_version], [19.1], [ ompimp="pgi" ], [ ompimp="intel" ], [ ompimp="intel" ] )
      		    mkl_lib="mkl_${ompimp}_lp64"
      		    mkl_omp="mkl_${ompimp}_thread"
      		    add_mkl_flag="-pgf90libs"
	       ;;
	    esac
            try_libdirs="$libdirs $MKLROOT/lib/intel64 $ld_library_path"
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
                    LDFLAGS="$add_mkl_flag $test_ldflags $try_loption"
		    # LIBS=""
                    # not sure the above is needed
                    if test "$use_openmp" -eq 0; then
		       # test MKL (no OMP)
 			      AC_SEARCH_LIBS(dgemm, $mkl_lib,
                             have_blas=1 have_mkl=1 
                             blas_libs="$try_loption $LIBS -lmkl_sequential -lmkl_core"
                             ldflags="$add_mkl_flag $ldflags",
                             echo "MKL not found",
                             -lmkl_sequential -lmkl_core $add_mkl_lib)
                    else
		       # test MKL (OMP)
 			      AC_SEARCH_LIBS(dgemm, $mkl_lib,
                             have_blas=1 have_mkl=1 
                             blas_libs="$try_loption $LIBS -l$mkl_omp -lmkl_core"
                             ldflags="$add_mkl_flag $ldflags",
                             echo "MKL not found",
                             -l$mkl_omp -lmkl_core $add_mkl_omp)
                    fi
                    if test "$ac_cv_search_dgemm" != "no"
                    then break ; fi
            done
            ;;
	    
    ppc64:* )
            #
            # search for ESSL - newer (?) powerPC machines
	    #
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

    ppc64-*:*  )
            #
            # assume ESSL without testing - old powerPC machines, BlueGene
	    #
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
	    
    arm:armflang )
	    # search for ARM libs - ARM compiler
            if test "$use_openmp" -eq 0; then 
               FFLAGS="-armpl"
            else 
               FFLAGS="-fopenmp -armpl=parallel" 
            fi 
            AC_SEARCH_LIBS(dgemm, armpl_arm,
                                   have_blas=1 have_armpl=1
                                   blas_libs=""
                                   ldflags="$ldflags \$(FFLAGS)",
                                   echo "armpl not found",
                                   yes)
            if test "$have_armpl" -eq 1; then 
               if test "$use_openmp" -eq 0; then 
                  fflags="$fflags  -armpl"
               else 
                  fflags="$fflags -armpl=parallel" 
               fi
            fi 
           ;;

    arm:gfortran )
	    # search for ARM libs - gfortran compiler
          try_libdirs="$libdirs $ARMPL_LIBRARIES $ld_library_path" 
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
                    LDFLAGS="$test_ldflags $try_loption"
		    # LIBS=""
                    # not sure the above is needed
                    #
                    if test "$use_openmp" -eq 0; then
 			      AC_SEARCH_LIBS(dgemm, armpl, 
                             have_blas=1 have_armpl=1 
                             blas_libs="$try_loption $LIBS "
                             ldflags="$ldflags",
                             echo "armpl not found",
                             )
			else
 			      AC_SEARCH_LIBS(dgemm, armpl_mp, 
                             have_blas=1 have_armpl=1 
                             blas_libs="$try_loption $LIBS "
                             ldflags="$ldflags",
                             echo "armpl  not found",
                             )
                    fi
                    if test "$ac_cv_search_dgemm" != "no"
                    then break ; fi
          done       
          ;;

    # obsolescent or obsolete architectures
    
    crayxt*:* )
            # check for acml - OBSOLETE?
            try_libdirs="$ld_library_path $libdirs"
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
            try_libdirs="$libdirs /usr/local/lib $ld_library_path"

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

    if test "$have_blas" -eq 0
    then
            # check for blas (in several directories)
            try_libdirs="$libdirs /usr/local/lib $ld_library_path"

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
fi

if test "$have_blas" -eq 0  ; then
    # No blas library found: use internal one (in lapack)
    blas_libs="\$(TOPDIR)/LAPACK/libblas.a" 
else
    echo setting BLAS_LIBS... $blas_libs
fi
blas_line="BLAS_LIBS=$blas_libs" 

AC_SUBST(blas_libs)
AC_SUBST(blas_line)
  
]
)
