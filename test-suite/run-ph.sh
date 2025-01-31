#!/bin/bash
#
# Copyright (C) 2024 Quantum ESPRESSO
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

# echo $0" "$@
if [[ "$1" == "1" ]]
then
  echo "Running PW ..."
# echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "11" ]]
then
  if [[ -e CRASH ]]
  then
    cat CRASH > $3
  else
    echo "Running PH ..."
    ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/ph.x ${PARA_SUFFIX} < $2 > $3 2> $4
    if [[ -e CRASH ]]
    then
      cat $3
    fi
  fi
elif [[ "$1" == "2" ]]
then
  echo "Running PH ..."
# echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/ph.x ${PARA_SUFFIX} < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/ph.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "3" ]]
then
  echo "Running Q2R ..."
# echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/q2r.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/q2r.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "4" ]]
then
  echo "Running MATDYN ..."
# echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/matdyn.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/matdyn.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "5" ]]
then
  echo "Running LAMBDA ..."
# echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/lambda.x < $2 > $3 2> $4"
  ${ESPRESSO_BUILD}/bin/lambda.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "6" ]]
then
  echo "Running DVSCF_Q2R ..."
# echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/dvscf_q2r.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/dvscf_q2r.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "7" ]]
then
  echo "Running POSTAHC ..."
# echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/postahc.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/postahc.x < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "8" ]]
then
  echo "Running MATDYN ..."
# echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/matdyn.x < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/matdyn.x < $2 > $3 2> $4
  cp matdyn.modes $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi
elif [[ "$1" == "9" ]]
then
   echo "Running DYNMAT ... "
   ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/dynmat.x < $2 > $3 2> $4
elif [[ "$1" == "12" ]]
then
  echo "Running PW ..."
  echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $3
  fi
  echo "Before running the multipole.py Python script, pip install numpy and spglib"
  echo " "
  pip install numpy
  pip3 install numpy
  pip install scipy
  pip3 install scipy
  pip install spglib
  pip3 install spglib
  echo " "
  echo "Running Python multipole.py preprocessing..."
  python3 multipole.py -e --order 3 --epsil_order 4 -p --mesh 3 3 3 --mesh_step 0.01 --alat 8.237B --ir_q > preprocessing.out
elif [[ "$1" == "13" ]]
then
  echo "Running PH ..."
  echo "${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/ph.x ${PARA_SUFFIX} < $2 > $3 2> $4"
  ${PARA_PREFIX} ${ESPRESSO_BUILD}/bin/ph.x ${PARA_SUFFIX} < $2 > $3 2> $4
  echo "Running Python postprocessing.."
  python3 multipole.py -f --order 3 --epsil_order 4 --alat 8.237B > postprocessing.out
  echo "postprocessing.out" >> $3
  cat postprocessing.out >> $3
  echo "epsilon.fmt" >> $3
  cat epsilon.fmt >> $3
  echo "born_charge.fmt" >> $3
  cat born_charge.fmt >> $3
  echo "quadrupole.fmt" >> $3
  cat quadrupole.fmt >> $3
  if [[ -e CRASH ]]
  then
    cat $3
  fi
fi

