# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_MPI], [

have_mpi=0
parallel=0

# some architectures require to link mpi libraries explicitly
F77=$mpif90 # use parallel compiler
if test "$mpi_libs" = ""
then
        # check directories in LD_LIBRARY_PATH too
        # (maybe they are already searched by default, but I'm not sure)
        ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`

        if test "$use_parallel" -ne 0
        then
                if test "$have_mpi" -eq 0
                        # check for mpi
                then
                        unset ac_cv_search_mpi_init # clear cached value
                        LDFLAGS="$test_ldflags"
                        LIBS="$mpi_libs"
                        AC_SEARCH_LIBS(mpi_init, mpi, 
                                       have_mpi=1 parallel=1 mpi_libs="$LIBS" try_dflags="$try_dflags -D__MPI")
                fi
        fi
else
        if test "$use_parallel" -ne 0
        then
                have_mpi=1
                parallel=1
                try_dflags="$try_dflags -D__MPI"
        fi
fi

# Configuring output message
if test "$mpi_libs" != "" ; then
   mpi_line="MPI_LIBS=$mpi_libs"
else
   mpi_line="@delete@"
fi

# Parallel report
if test "$use_parallel" -ne 0
then
        if test "$parallel" -ne 0
        then
                parallel_report="Parallel environment detected successfully.\\
Configured for compilation of parallel executables."
        else
                parallel_report="Parallel environment not detected \
\(is this a parallel machine?\).\\
Configured for compilation of serial executables."
        fi
else
        parallel_report="Configured for compilation of serial executables."
fi

  AC_SUBST(mpi_libs)
  AC_SUBST(mpi_line)
  AC_SUBST(parallel_report)
  
  ]
)
