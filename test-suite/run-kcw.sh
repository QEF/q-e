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

#echo $0" "$@
# First SCF, no k pools
if [[ "$1" == "1" ]]
then
  echo "Running PW ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi
#rm -f input_tmp.in
if [[ "$1" == "2" ]]
then
  echo "Running Wannier90 PP ..."
  #echo "${ESPRESSO_ROOT}/bin/wannier90.x -pp $2 "
  ${ESPRESSO_ROOT}/bin/wannier90.x -pp $2 
  cp *nnkp $3
  if [[ -e Si.werr ]]
  then
    cat Si.werr
  fi
fi
if [[ "$1" == "3" ]]
then
  echo "Running pw2wannier90 ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw2wannier90.x < $2 > $3 2>$4 "
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw2wannier90.x < $2 > $3 2>$4
  cp *eig $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi
if [[ "$1" == "4" ]]
then
  echo "Running Wannier90 ..."
  #echo "${ESPRESSO_ROOT}/bin/wannier90.x $2 "
  ${ESPRESSO_ROOT}/bin/wannier90.x  $2 
  cp *.wout $3
  if [[ -e Si.werr ]]
  then
    cat Si.werr
  fi
fi

if [[ "$1" == "5" ]]
then
  echo "Running W90-KCW interface ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcw.x ${PARA_SUFFIX}< $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

if [[ "$1" == "15" ]]
then
  echo "Running W90-KCW interface ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcw.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcw.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

if [[ "$1" == "6" ]]
then
  echo "Running KCW screen ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcw.x ${PARA_SUFFIX} < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

if [[ "$1" == "16" ]]
then
  echo "Running KCW screen ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcw.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcw.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

if [[ "$1" == "7" ]]
then
  echo "Running KCW hamiltonian ..."
  #echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcw.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcw.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi
