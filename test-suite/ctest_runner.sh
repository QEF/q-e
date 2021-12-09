#!/bin/bash

# helper script to facilitate ctest
if [ $# -lt 3 ] ; then
  echo "Need at least 3 arguments"
  echo "ctest_runner.sh label outfile exe [args...]"
  echo "  label   : enter pw for pw.x. The label is converted to upper cases and used for read environment variable, for example CTEST_PW_OPTIONS"
  echo "  outfile : name of the file which repeats the standard output"
  echo "  exe     : excutable, full path preferred. It may not be a qe executable but an MPI runner"
  exit 1
fi

if [ ! -z "$CTEST_RESOURCE_GROUP_0_NVIDIA_GPUS" ] ; then
  GPU_ID=`echo $CTEST_RESOURCE_GROUP_0_NVIDIA_GPUS | sed "s/id:\(.*\),.*$/\1/"`
  echo "Assign GPU $GPU_ID to the run"
  export CUDA_VISIBLE_DEVICES=$GPU_ID
fi

ARGS=""

for i in `seq 3 $#`
do
  eval "arg=\${$i}"
  ARGS=$ARGS" "$arg
done

# access additional arguments from environment variables ${LABEL}_OPTIONS
label=$1
LABEL=`echo ${label%.*} | tr '[:lower:]' '[:upper:]'`
#echo LABEL $LABEL
option_var=CTEST_${LABEL}_OPTIONS

if [ ! -z "${!option_var}" ] ; then
  ARGS=$ARGS" "${!option_var}
fi

echo "Full test run command is $ARGS"
$ARGS | tee $2
exit ${PIPESTATUS[0]}
