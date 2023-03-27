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
  export PARA_SUFFIX="-npool $QE_USE_MPI"
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

echo $0" "$@
if [[ "$1" == "0" ]]
then
  echo "Running PW ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "1" ]]
then
  echo "Running PW ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "2" ]]
then
  echo "Running PH ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/ph.x ${PARA_SUFFIX} -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/ph.x ${PARA_SUFFIX} -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  echo "Gather results in save"
  python3 ../../EPW/bin/pp.py < pp.in
elif [[ "$1" == "3" ]]
then
  echo "Running EPW ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/epw.x ${PARA_SUFFIX} -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/epw.x ${PARA_SUFFIX} -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "4" ]]
then
  echo "Running Q2R ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/q2r.x -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/q2r.x -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  echo "Gather results in save"
  python3 ../../EPW/bin/pp.py < pp.in
elif [[ "$1" == "5" ]]
then
  echo "Removing restart files ..."
  echo "Running EPW ..."
######  rm *.Fin_restart1 *.Fin_restartcb1 restart.fmt
  rm restart.fmt
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/epw.x ${PARA_SUFFIX} -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/epw.x ${PARA_SUFFIX} -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

#rm -f input_tmp.in
