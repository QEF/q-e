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
echo $1
if [[ "$1" == "0" ]]
then
  echo "Running PW ..."
  echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "1" ]]
then
  echo "Running PW ..."
  echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "2" ]]
then
  echo "Running PH ..."
  echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/ph.x ${PARA_SUFFIX} -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/ph.x ${PARA_SUFFIX} -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  echo "Gather results in save"
  python3 ../../EPW/bin/pp.py < pp.in
elif [[ "$1" == "3" ]]
then
  echo "Running EPW ..."
  echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/epw.x ${PARA_SUFFIX} -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/epw.x ${PARA_SUFFIX} -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "4" ]]
then
  echo "Running Q2R ..."
  echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/q2r.x -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/q2r.x -input $2 > $3 2> $4
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
  echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/epw.x ${PARA_SUFFIX} -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/epw.x ${PARA_SUFFIX} -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "6" ]]
then
  echo "Running MATDYN ..."
# echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/matdyn.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/matdyn.x < $2 > $3 2> $4
  cp matdyn.modes $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "7" ]]
then
  # ph.x, electron_phonon = 'ahc'
  rm -rf save/ahc_dir/
  echo "Running PH ..."
  echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/ph.x ${PARA_SUFFIX} -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/ph.x ${PARA_SUFFIX} -input $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "8" ]]
then
  echo "Running POSTAHC ..."
# echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/postahc.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/postahc.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "9" ]]
then
  # nscf2supercond.x
  echo "Running NSCF2SUPERCOND ..."
  echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/nscf2supercond.x ${PARA_SUFFIX} -input $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/nscf2supercond.x ${PARA_SUFFIX} -input $2 > $3 2> $4
  cat *.bands.*.dat >> $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

#rm -f input_tmp.in
