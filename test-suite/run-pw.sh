#!/bin/bash
#
# Copyright (C) 2025 Quantum ESPRESSO
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

if [[ "$QE_USE_MPI" != "" ]]; then
  export PARA_PREFIX="mpirun -np $QE_USE_MPI"
  export PARA_SUFFIX=" "
  if [[ "$QE_USE_BGRP" != "" ]]; then
    export PARA_SUFFIX=" -nb $QE_USE_BGRP "
  fi
elif [[ "$QE_USE_BGRP" != "" ]]; then
  echo "WARNING: ignoring NBGRP because NPROCS is not set"
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi
if [[ "$1" == "1" ]]
then
  #echo "Running PW ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "2" ]]
then
  if [[ -e CRASH ]]
  then 
    cat CRASH > $3
  else
    #echo "Running PW ..."
    #echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4"
    ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
    if [[ -e CRASH ]]
    then
      cat $3
    fi
  fi 
elif [[ $1 == "48" ]]
then
  export PARA_SUFFIX="$PARA_SUFFIX --pw2casino"
  # echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} -input $1 > $2 2> $3"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} -input $2 > $3 2> $4
elif [[ $1 == "22" ]]
then
  # This is a restart test, need to clean up previous results if present
  rm -rf md_restart_verlet.save
  cp md_restart_verlet_original.md md_restart_verlet.md
  # echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} -input $1 > $2 2> $3"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} -input $2 > $3 2> $4
else
  # no arg provided input is $1 output $2 stderr is $3 
  # echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} -input $1 > $2 2> $3"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} -input $1 > $2 2> $3
fi

rm -f input_tmp.in
