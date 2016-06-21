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

#include ${ESPRESSO_ROOT}/test-suite/ENVIRONMEN
bash ../ENVIRONMENT

if [[ $QE_USE_MPI == 1 ]]; then
  export PARA_PREFIX="mpirun -np ${TESTCODE_NPROCS}"
else
  unset PARA_PREFIX
fi

echo $0" "$@
if [[ "$1" == "1" ]]
then
  echo "Running PW ..."
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x < $2 > $3 2> $4
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x < $2 > $3 2> $4"
elif [[ "$1" == "2" ]]
then
  echo "Running PH ..."
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/ph.x < $2 > $3 2> $4
  echo "gather results in save" 
  python pp.py < pp.in
elif [[ "$1" == "3" ]]
then
  echo "Running EPW ..."
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/epw.x < $2 > $3 2> $4
fi

#rm -f input_tmp.in
