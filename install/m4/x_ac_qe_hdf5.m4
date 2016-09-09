# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_HDF5], [

AC_ARG_WITH(hdf5,
   [AS_HELP_STRING([--with-hdf5],
       [(no|<path>) Use HDF5, a valid <path> must be specified (default: no)])],
   [if  test "$withval" = "no" ; then
      with_hdf5=0
   else
      with_hdf5_path="$withval"
      with_hdf5=1
   fi],
   [with_hdf5=0])

hdf5_libs=""
have_hdf5=0

if test "$use_parallel" -ne 0; then

  if test "$with_hdf5" -ne 0 && test "$with_hdf5_path" != "yes"; then

    # Checking compiler compatibility: GCC >= 4.9
    if test "x$f90_in_mpif90" = xgfortran && 
       test "$f90_major_version" -le "4" && 
       test "$f90_minor_version" -lt "9" ; then

      AC_MSG_RESULT(no)
      AC_MSG_WARN([*** HDF5 support requires GNU GFORTRAN >= 4.9 ])

    else

      # Test if it is really installed where it has been specified
      AC_LANG_POP(Fortran 77)
      AC_LANG_PUSH(C)

      try_libdirs="$with_hdf5_path/lib"
      for dir in $try_libdirs
      do
          unset ac_cv_search_H5Fcreate

          if test "$dir" = "none"
          then
            try_loption=
          else
            try_loption="-L$dir"
          fi

          FFLAGS="$test_fflags"
          LDFLAGS="$test_ldflags $try_loption"
          #LIBS="-lhdf5"
          LIBS=""

          AC_SEARCH_LIBS(H5Fcreate, hdf5_fortran, [have_hdf5=1])

          if test "$ac_cv_search_H5Fcreate" != "no"
          then break ; fi
      done

      AC_LANG_POP(C)
      AC_LANG_PUSH(Fortran 77)

      if test "$have_hdf5" -eq 1 ; then
          AC_CHECK_FILE($with_hdf5_path/include/hdf5.mod,,[
              AC_MSG_WARN([***HDF5 Fortran extensions not found])
              have_hdf5=0])
      fi
      
      if test "$have_hdf5" -eq 1 ; then
          hdf5_libs="-L$with_hdf5_path/lib -lhdf5_fortran -lhdf5"
          try_iflags="$try_iflags -I$with_hdf5_path/include"
          try_dflags="$try_dflags -D__HDF5"
      fi

      hdf5_line="HDF5_LIBS=$hdf5_libs"
    fi
  fi
#else
#    AC_MSG_WARN([HDF5 support is for parallel execution only])
fi

  AC_SUBST(hdf5_libs)
  AC_SUBST(hdf5_line)
  ]
)
