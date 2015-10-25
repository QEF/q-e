# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_MASS], [

  AC_MSG_CHECKING([MASS])
 
 # check for mass on aix
if test "$mass_libs" = ""
then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        case "$arch" in
        aix | ppc64-bg )
                # check for mass (in several directories)
                try_libdirs="/opt/ibmcmp/xlmass/bg/7.3/bglib64 /opt/ibmcmp/xlmass/bg/4.4/bglib /cineca/lib /cineca/lib/mass"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

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
                        LIBS=""
                        AC_SEARCH_LIBS(vexp, massvp4 massv, , , -lmass)
                        if test "$ac_cv_search_vexp" = "-lmassvp4" \
                                -o "$ac_cv_search_vexp" = "-lmassv"
                        then mass_libs="$try_loption $ac_cv_search_vexp -lmass"
                        fi
                        if test "$ac_cv_search_vexp" != "no" ; then break ; fi
                done
                ;;
        ppc64-bgq )
                # check for mass (in several directories)
                try_libdirs="/opt/ibmcmp/xlmass/bg/7.3/bglib64"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

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
                        LIBS=""
                        AC_SEARCH_LIBS(vexp, massv, , , -lmass_simd)
                        if test "$ac_cv_search_vexp" = "-lmassv"
                        then mass_libs="$try_loption $ac_cv_search_vexp -lmass_simd"
                        fi
                        if test "$ac_cv_search_vexp" != "no" ; then break ; fi
                done
                ;;

        ppc64* )
                # check for mass (in several directories)
                try_libdirs="/usr/local/lib /opt/ibmcmp/xlmass/*/lib64"
                try_libdirs="$libdirs $try_libdirs $ld_library_path"

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
                        LIBS=""
                        AC_SEARCH_LIBS(vexp, massvp4_64, , , -lmass_64)
                        if test "$ac_cv_search_vexp" = "-lmassvp4_64"
                        then mass_libs="$try_loption $ac_cv_search_vexp -lmass_64"
                        fi
                        if test "$ac_cv_search_vexp" != "no" ; then break ; fi
                done
                ;;

        esac
fi

if test "$mass_libs" != ""; then 
   try_dflags="$try_dflags -D__MASS"
   if test "$arch" = "ppc64-bg"; then
   # BlueGene wants this when mass libs are loaded, SP6 doesn't want this!
     ldflags="$ldflags -Wl,--allow-multiple-definition"
   fi
   if test "$arch" = "ppc64-bgq"; then
   # BlueGene wants this when mass libs are loaded, SP6 doesn't want this!
     ldflags="$ldflags -Wl,--allow-multiple-definition"
   fi
fi

# Configuring output message
if test "$mass_libs" != "" ; then
   mass_line="MASS_LIBS=$mass_libs"
else
   mass_line="@delete@"
fi

  AC_MSG_RESULT(${mass_libs})
  
  AC_SUBST(mass_libs)
  AC_SUBST(mass_line)
  
  ]
)
