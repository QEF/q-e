#!/bin/bash
#
# Copyright (C) 2020 Quantum ESPRESSO Foundation
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

if [[ "$QE_USE_MPI" != "" ]]; then
  export PARA_PREFIX="mpirun -np $QE_USE_MPI"
  export PARA_SUFFIX=" "
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

if [[ "$1" == "1" ]]
then
  # echo "Running PW ..."
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi  
elif [[ "$1" == "2" ]]
then
  # echo "Running PP ..."
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/ppacf.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

