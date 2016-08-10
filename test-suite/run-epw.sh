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
  export PARA_SUFFIX="-npool ${TESTCODE_NPROCS}"
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

echo $0" "$@
if [[ "$1" == "1" ]]
then
  echo "Running PW ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi  
elif [[ "$1" == "2" ]]
then
  echo "Running PH ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/ph.x < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/ph.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  echo "Gather results in save" 
  python pp.py < pp.in
elif [[ "$1" == "3" ]]
then
  echo "Running EPW ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/epw.x ${PARA_SUFFIX} < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/epw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

#rm -f input_tmp.in
