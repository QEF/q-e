#!/bin/bash


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
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} < $2 > $3 2> $4
  if [[ -e CRASH ]]
  then
    cat $4
  fi
elif [[ "$1" == "2" ]]
then
  echo "Running PW OSCDFT ..."
  cp "$2" oscdft.in
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} -oscdft < $3 > $4 2> $5
  if [[ -e CRASH ]]
  then
    cat $5
  fi
elif [[ "$1" == "3" ]]
then
  echo "Running OSCDFT_PP ..."
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/pw.x ${PARA_SUFFIX} < $3 > "$3".out.tmp 2> "$3".err.tmp
  cp "$2" oscdft.in
  ${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/oscdft_pp.x ${PARA_SUFFIX} < $4 > $5 2> $6
  if [[ -e CRASH ]]
  then
    cat $6
  fi
fi

rm -f input_tmp.in
