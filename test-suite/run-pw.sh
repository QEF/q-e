#!/bin/bash
#
# Copyright (C) 2001-2016 Quantum ESPRESSO group
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Maintainer: Filippo Spiga (filippo.spiga@quantum-espresso.org)

source ${ESPRESSO_ROOT}/test-suite/ENVIRONMENT


# Additional stuff before run special test-cases

if test "$2" = "vdw1.in" || test "$1" = "vdw2.in" ; then
   if ! test -f ${ESPRESSO_PSEUDO}/vdW_kernel_table ; then
      echo -n "Generating kernel table - May take several minutes..."
      ${ESPRESSO_ROOT}/PW/src/generate_vdW_kernel_table.x 
      mv vdW_kernel_table ${ESPRESSO_PSEUDO}/
      echo "kernel table generated in ${ESPRESSO_PSEUDO}/vdW_kernel_table"
   fi
fi

if test "$2" = "vdw6.in" ; then
   if ! test -f ${ESPRESSO_PSEUDO}/rVV10_kernel_table ; then
      echo -n "Generating kernel table - May take several minutes..."
      ${ESPRESSO_ROOT}/PW/src/generate_rVV10_kernel_table.x
      mv rVV10_kernel_table ${ESPRESSO_PSEUDO}/
      echo "kernel table generated in ${ESPRESSO_PSEUDO}/rVV10_kernel_table"
   fi
fi

# -- Uncomment to run in parallel
#export PARA_PREFIX="mpirun -np 4"
unset PARA_PREFIX

${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x "$@"

rm -f input_tmp.in
