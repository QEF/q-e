#!/bin/bash
#
#For serial
set -x


echo "Running PIOUD (serial)... $@"
${ESPRESSO_ROOT}/bin/pioud.x ${PARA_SUFFIX} -in  $1 > $2 2> $3
rm -f pw_?.in

if grep -Eq "nbeadMD *= *[2-9]" $1; then
  NIMAGE=2
else
  NIMAGE=1
fi

#For parallel
if [[ $QE_USE_MPI == 1 ]]; then
  export PARA_PREFIX="mpirun -np ${TESTCODE_NPROCS}"
  export PARA_SUFFIX="-nimage $NIMAGE"
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

echo "Running PIOUD (parallel)..."
${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pioud.x ${PARA_SUFFIX} -in $1 > $2 2> $3
rm -f pw_?.in