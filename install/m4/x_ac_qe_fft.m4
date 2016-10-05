# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_FFT], [

have_fft=0
have_fft_include=0

  AC_MSG_CHECKING([FFT])
 
# check for FFT libraries (no check for explicit openmp)
# supported vendor replacements:
#   essl on aix and some IBM linux machines
#   SUNperf on sparc
#   ASL/Mathkeisan on Nec
#   acml on amd
if test "$fft_libs" = "" && test "$use_openmp" -eq 0
then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        case "$arch" in
        aix )
                # check for essl
                unset ac_cv_search_dcft # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$fft_libs"
                AC_SEARCH_LIBS(dcft, essl, have_fft=1 fft_libs="$LIBS")
            ;;
        ppc64 | ppc64-mn )
                # check for essl
                unset ac_cv_search_dcft # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$fft_libs"
                AC_SEARCH_LIBS(dcft, essl, have_fft=1 fft_libs="$LIBS")
            ;;
        ppc64-bg | ppc64-bgq )
                # check for esslbg
                unset ac_cv_search_dcft # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$fft_libs $blas_libs"
                AC_SEARCH_LIBS(dcft, esslbg, have_fft=1 fft_libs="$LIBS")
            ;;
        sparc )
                # check for SUNperf FFT library on Sun Sparcs
                # but not on solaris PC! it is slower than FFTW
                unset ac_cv_search_zfft3i # clear cached value
                FFLAGS="$test_fflags"
                LDFLAGS="$test_ldflags"
                LIBS="$libs"
                AC_SEARCH_LIBS(zfft3i, sunperf, have_fft=1
                               try_dflags="$try_dflags -D__SUNPERF"
                               fft_libs="-xlic_lib=sunperf $LIBS")
                ;;
        necsx )
                if test "$use_fft_mathkeisan" -ne 0
                then
                   #sx5-nec or sx6-nec or sx8-nec: check in (/SX)/usr/lib
                   #sx8-nec-idris: check in /SX/opt/mathkeisan/inst/lib0
                   try_libdirs="/SX/usr/lib /SX/opt/mathkeisan/inst/lib0"
                   #check for Mathkeisan (Cray simple precision )
                   #search for initialization subroutine
                   echo $ECHO_N  "Searching in Mathkeisan" $ECHO_C
                   for dir in none $try_libdirs
                   do
                        unset ac_cv_search_zftfax # clear cached value
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
                        AC_SEARCH_LIBS(zftfax, fft, have_fft=1
                             try_dflags="$try_dflags try_dflags_fft_mathkeisan"
                                       fft_libs="$try_loption $LIBS")
                        if test "$ac_cv_search_zftfax" != "no"
                        then break ; fi
                   done
                fi
                if test "$use_fft_asl" -ne 0
                then
                   #check for asl in (/SX)/usr/lib
                   try_libdirs="/SX/usr/lib"
                   #search for initialization subroutine
                   echo $ECHO_N  "Searching in Asl" $ECHO_C
                   for dir in none $try_libdirs
                   do
                        unset ac_cv_search_zfc3cl # clear cached value
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
                        AC_SEARCH_LIBS(zfc3cl, asl, have_fft=1
                             asl_libs="$try_loption $LIBS"
                             try_dflags="$try_dflags $try_dflags_fft_asl"
                             fft_libs="$fft_libs $asl_libs")
                        if test "$ac_cv_search_zfc3cl" != "no"
                        then break ; fi
                   done
                fi
                if test "$use_fft_para" -ne 0
                then
                   try_dflags="$try_dflags $try_dflags_fft_para"
                fi
                ;;
        esac

fi

# Check for fftw v3, both native and MKL interfaces
if test "$have_fft" -eq 0 
then

        try_libdirs="/usr/local/lib"
        try_libdirs="$libdirs $try_libdirs $ld_library_path "

        AC_LANG_POP(Fortran 77)
        AC_LANG_PUSH(C)
        #for dir in none $try_libdirs
        #do
        unset ac_cv_lib_mkl_intel_lp64_DftiComputeForward
        if test "$dir" = "none"
        then
                try_loption=
        else
                echo $ECHO_N "in $dir: " $ECHO_C
                try_loption="-L$dir"
        fi

        CFLAGS="$test_cflags"
        CPPFLAGS="$test_cppflags"
        LDFLAGS=" $test_ldflags $try_loption"
        LIBS="$fft_libs"

        #here we check if dfti explicit calls work
        #it should work with blas flags
        AC_CHECK_LIB([mkl_intel_lp64],DftiComputeForward,have_fft=1,
        ,$blas_libs -lm)

        if test "$have_fft" == "1" 
        then
              try_incdir=""
              try_incdir="$try_incdir /opt/intel/Compiler/*/*/mkl/include
                          /opt/intel/mkl/*/include
                          /opt/intel/mkl*/include"
              try_incdir="$try_incdir $MKLROOT/include $MKL_INCLUDE $CPATH $FPATH"
              try_incdir="`echo $try_incdir | sed 's/:/ /g'`"
              # 
              for inc in $try_incdir
              do
                 if test "$cross_compiling" == "no"
                 then
                    # when cross-compilation is disabled I can check file paths... 
                    AC_CHECK_FILE( $inc/mkl_dfti.f90, 
                                have_fft_include=1, have_fft_include=0) 
                 else
                    # when cross-compilation is enabled I need to rely on header 
                    # checking (bit it complains...
                    AC_CHECK_HEADER($inc/mkl_dfti.f90,
                                have_fft_include=1,have_fft_include=0)
                 fi
                 if test "$have_fft_include" == "1"
                 then
                   try_iflags="$try_iflags -I$inc"
                   try_dflags="$try_dflags -D__DFTI"
                   have_fft=1
                   break
                 else
                   # continue the search
                   have_fft=0
                 fi
              done
        fi

        AC_LANG_POP(C)
        AC_LANG_PUSH(Fortran 77)

        if test "$have_fft" -eq 0
        then
          for dir in none $try_libdirs
          do
                  unset ac_cv_search_dfftw_execute_dft # clear cached value
                  if test "$dir" = "none"
                  then
                        try_loption=
                  else
                        echo $ECHO_N "in $dir: " $ECHO_C
                        try_loption="-L$dir"
                  fi

                  CFLAGS="$test_cflags"
                  CPPFLAGS="$test_cppflags"
                  LDFLAGS=" $test_ldflags $try_loption"
                  LIBS="$fft_libs"

                  if test "$use_openmp" -eq 1
                  then
                    AC_SEARCH_LIBS(dfftw_execute_dft, fftw3_omp, have_fft=1
                               fft_libs="$try_loption $LIBS -lfftw3", , -lfftw3 -lm)
                  else
                    AC_SEARCH_LIBS(dfftw_execute_dft, fftw3, have_fft=1
                               fft_libs="$try_loption $LIBS", , -lm)
                  fi

                  if test "$have_fft" -eq 1
                  then
                        try_dflags="$try_dflags -D__FFTW3"
                        try_incdir="$FFTW_INCLUDE $FFTW_INC $INCLUDE_PATH $CPATH $FPATH"
                        for inc in $try_incdir
                        do
                           #AC_LANG_POP([Fortran 77])
                           #AC_LANG_PUSH([C])

                           AC_COMPILE_IFELSE([include "fftw3.f03"],have_fft_include=1,)
                           
                           #AC_LANG_POP([C])
                           #AC_LANG_PUSH([Fortran 77])

                           if test "$have_fft_include" -eq 1
                           then
                             try_iflags="$try_iflags -I$inc"
                             break
                           fi
                        done
                        break
                  fi

          done
        fi
fi

  AC_MSG_RESULT(${fft_libs})

fft_line="FFT_LIBS=$fft_libs"

# if no valid FFT library was found, use the local copy
# (This happens also if OpenMP is enabled...)
if test "$have_fft" -eq 0
then
        case "$arch" in
        ppc64-bg | ppc64-bgq )
            try_dflags="$try_dflags -D__LINUX_ESSL"
            ;;
        * )
		    try_dflags="$try_dflags -D__FFTW"
            ;;
        esac
fi

AC_SUBST(fft_libs)
AC_SUBST(fft_line)
  
  ]
)
