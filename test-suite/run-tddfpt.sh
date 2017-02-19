#!/bin/bash
#
# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

bash ../ENVIRONMENT

if [[ $QE_USE_MPI == 1 ]]; then
  export PARA_PREFIX="mpirun -np ${TESTCODE_NPROCS}"
  export PARA_SUFFIX=" "
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
  echo "Running TURBO LANCZOS ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/turbo_lanczos.x < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/turbo_lanczos.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "3" ]]
then
  echo "Running TURBO SPECTRUM CHI ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/turbo_spectrum.x ${PARA_SUFFIX} < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/turbo_spectrum.x ${PARA_SUFFIX} < $2 > $3 2> $4
  cp *.plot_chi.dat $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "4" ]]
then
  echo "Running TURBO SPECTRUM EELS ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/turbo_spectrum.x ${PARA_SUFFIX} < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/turbo_spectrum.x ${PARA_SUFFIX} < $2 > $3 2> $4
  cp *.plot_eps.dat $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "5" ]]
then
  echo "Running TURBO EELS ..."
  echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/turbo_eels.x ${PARA_SUFFIX} < $2 > $3 2> $4"  
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/turbo_eels.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

