#!/bin/bash
#
# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

if [[ "$QE_USE_MPI" != "" ]]; then
  export PARA_PREFIX="mpirun -np $QE_USE_MPI"
else
  unset PARA_PREFIX
fi

echo $0" "$@
if [[ "$1" == "0" ]]
 then
 echo "Running ZG ..."
 echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/ZG.x -input $2 > $3 2> $4"
 ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/ZG.x -input $2 > $3 2> $4
fi
#
if [[ -e CRASH ]]
then
    cat $2
fi  
