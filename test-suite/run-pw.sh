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
  export PARA_SUFFIX=" "
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

if [[ "$1" == "1" ]]
then
  echo "Running PW ..."
# echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "2" ]]
then
  if [[-e CRASH ]]
  then 
    cat CRASH > $3
  else
    echo "Running PW ..."
# echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4"
    ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
    if [[ -e CRASH ]]
    then
      cat $3
    fi
  fi 
elif [[ "$1" = "plugin-pw2casino_1.in" ]] || [[ "$1" = "plugin-pw2casino_2.in" ]]
then
  export PARA_SUFFIX="$PARA_SUFFIX --pw2casino"
  # echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} -input $1 > $2 2> $3"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} -input $1 > $2 2> $3
elif [[ "$1" = "md_restart_verlet.in" ]]
then
  # This is a restart test, need to clean up previous results if present
  rm -rf md_restart_verlet.save
  cp md_restart_verlet_original.md md_restart_verlet.md
  # echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} -input $1 > $2 2> $3"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} -input $1 > $2 2> $3
else
  # echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} -input $1 > $2 2> $3"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} -input $1 > $2 2> $3
fi

rm -f input_tmp.in
