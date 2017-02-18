#!/bin/bash
#
# Copyright (C) 2001 Quantum ESPRESSO
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

include ${ESPRESSO_ROOT}/test-suite/ENVIRONMENT

if [ $QE_USE_MPI == 1 ]; then
  export PARA_PREFIX="mpirun -np ${TESTCODE_NPROCS}"
else
  unset PARA_PREFIX
fi

# Additional stuff before run special test-cases

if test "$2" = "vdw1.in" || test "$1" = "vdw2.in" ; then
   if ! test -f ${ESPRESSO_PSEUDO}/vdW_kernel_table ; then
      echo -n "Generating kernel table - May take several minutes..."
      ${PARA_PREFIX} ${ESPRESSO_ROOT}/PW/src/generate_vdW_kernel_table.x 
      mv vdW_kernel_table ${ESPRESSO_PSEUDO}/
      echo "kernel table generated in ${ESPRESSO_PSEUDO}/vdW_kernel_table"
   fi
fi

if test "$2" = "vdw6.in" ; then
   if ! test -f ${ESPRESSO_PSEUDO}/rVV10_kernel_table ; then
      echo -n "Generating kernel table - May take several minutes..."
      ${PARA_PREFIX} ${ESPRESSO_ROOT}/PW/src/generate_rVV10_kernel_table.x
      mv rVV10_kernel_table ${ESPRESSO_PSEUDO}/
      echo "kernel table generated in ${ESPRESSO_PSEUDO}/rVV10_kernel_table"
   fi
fi

${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x "$@"

rm -f input_tmp.in
