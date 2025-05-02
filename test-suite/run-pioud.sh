#!/bin/bash
#
#For serial
# set -x


# echo "Running PIOUD (serial)... $@"
# ${ESPRESSO_ROOT}/bin/pioud.x ${PARA_SUFFIX} -in  $1 > $2 2> $3
# rm -f pw_?.in

NIMAGE=1

#For parallel
if [[ "$QE_USE_MPI" != "" ]]; then
  if grep -Eq "nbeadMD *= *[2-9]" $1; then
    NIMAGE=2
  fi
  export PARA_PREFIX="mpirun -np $QE_USE_MPI"
  export PARA_SUFFIX=" -nimage $NIMAGE"
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

echo "Running PIOUD..."
${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pioud.x ${PARA_SUFFIX} -in $1 > $2 2> $3
rm -f pw_?.in
