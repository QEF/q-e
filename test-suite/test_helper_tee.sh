#!/bin/bash

if [ $# -lt 2 ] ; then
  echo "Need at least 2 arguments"
  exit 1
fi

ARGS=""

for i in `seq 2 $#`
do
  eval "arg=\${$i}"
  ARGS=$ARGS" "$arg
done

$ARGS | tee $1
