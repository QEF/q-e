# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_F90RULE], [

AC_PROG_MAKE_SET
echo $ECHO_N "checking whether Fortran files must be preprocessed... $ECHO_C"
if test "$have_cpp" -ne 0
then
        f90rule="\$(MPIF90) \$(F90FLAGS) -c \$<"
        echo "${ECHO_T}no"
else
        f90rule="\$(CPP) \$(CPPFLAGS) \$< -o \$(*)_tmp.f90 ; \\
        \$(MPIF90) \$(F90FLAGS) -c \$(*)_tmp.f90 -o \$(*).o"
        echo "${ECHO_T}yes"
fi

AC_SUBST(f90rule)

])
