# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_HDF5], [

  AC_MSG_CHECKING([hdf5])
 
AC_ARG_WITH(hdf5,
   [AS_HELP_STRING([--with-hdf5],
       [(yes|no|<path>) Use HDF5. Self-compile or a <path> can be specified (default: no)])],
   [if  test "$withval" = "yes" ; then
      with_hdf5=1
   elif  test "$withval" = "no" ; then
      with_hdf5=0
   else
      with_hdf5=2
      with_hdf5_path="$withval"
   fi],
   [with_hdf5=0])

hdf5_libs=""
hdf5_inc=""
have_hdf5=0
hdf5_libs_switch="disabled"

if test "$with_hdf5" -eq 2 ; then
     hdf5_libs="-L$with_hdf5_path/lib -lhdf5_fortran -lhdf5"
     try_iflags="$try_iflags -I$with_hdf5_path/include"
     try_dflags="$try_dflags -D__HDF5"
fi

if test "$with_hdf5" -eq 1 ; then
    try_libdirs="/usr/lib64"
    try_libdirs="$libdirs $try_libdirs $ld_library_path"

    have_hdf5=0

    AC_LANG_POP(Fortran 77)
    AC_LANG_PUSH(C)

    for dir in none $try_libdirs
    do
        # unset ac_cv_search_H5Fcreate ac_cv_search_hdfh 
        
        if test "$dir" = "none"
        then
            try_loption=
        else
            echo $ECHO_N "in $dir: " $ECHO_C
            try_loption="-L$dir"
        fi
        
        FFLAGS="$test_fflags"
        LDFLAGS="$test_ldflags $try_loption"
        LIBS="-lhdf5"

        # Detect the location Fortran module "USE HDF5" is the complex part 
        # AC_CHECK_HEADER(...)

        AC_SEARCH_LIBS(H5Fcreate, hdf5_fortran, 
          [hdf5_libs="$try_loption -lhdf5_fortran -lhdf5" have_hdf5=1])

        if test "$ac_cv_search_H5Fcreate" != "no"
        then break ; fi
    done

    AC_LANG_POP(C)
    AC_LANG_PUSH(Fortran 77)

    if test "$have_hdf5" -eq 1 ; then
        try_iflags="$try_iflags -I$with_hdf5_path/include"
        try_dflags="$try_dflags -D__HDF5"   
    fi
fi

  hdf5_line="HDF5_LIBS=$hdf5_libs"

  AC_SUBST(hdf5_libs)
  AC_SUBST(hdf5_line)
  ]
)
