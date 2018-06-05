#!/bin/bash
#
# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

if [ $QE_USE_MPI == 1 ]; then
  export PARA_PREFIX="mpirun -np ${TESTCODE_NPROCS}"
#
# available flags: 
#                  -ni n        number of images        (or -nimage)
#                               (only for NEB; for PHonon, see below)
#                  -nk n        number of pools         (or -npool, -npools)
#                  -nb n        number of band groups   (or -nbgrp,-nband_group)
#                  -nt n        number of task groups   (or -ntg, -ntask_groups)
#                  -nd n        number of processors for linear algebra 
#                                            (or -ndiag, -northo) 
#
  export PARA_POSTFIX=" -nk 1 -nd 1 -nb 1 -nt 1 "
else
  unset PARA_PREFIX
  unset PARA_POSTFIX
fi

echo ' RUNNING ',${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/cp.x ${PARA_POSTFIX} "$@"
${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/cp.x ${PARA_POSTFIX} "$@"

rm -f input_tmp.in
