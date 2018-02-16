#!/bin/bash

MPIEXEC="mpirun -np 3 "
GPU_ARCH=35
CUDA_RUNTIME=8.0

declare -a flags=()

add_cuda=""
add_cudampi=""

while getopts ":smcnh" opt; do
  case $opt in
    s)
      echo "Serial build scheduled!" >&2
      flags+=("")
      ;;
    m)
      echo "MPI build scheduled!" >&2
      flags+=(" -D__MPI")
      ;;
    c)
      echo "CUDA build scheduled!" >&2
      add_cuda="yes"
      ;;
    n)
      echo "CUDA+MPI build scheduled!" >&2
      add_cudampi="yes"
      ;;
    h)
      echo "-s : serial build" >&2
      echo "-m : mpi build" >&2
      echo "-c : cuda build" >&2
      echo "-n : gpu+mpi (nvlink) build" >&2
      echo "-h : this help" >&2
      echo "" >&2
      echo "Note: always use -c and -n with -s and/or -m" >&2
      exit
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit
      ;;
  esac
done
if [ "$add_cuda" != "" ]; then
  for flag in "${flags[@]}"
  do
    flags+=(" $flag -Mcuda=cc${GPU_ARCH},cuda${CUDA_RUNTIME} -D__CUDA")
  done
fi
if [ "$add_cudampi" != "" ]; then
  flags+=(" -Mcuda=cc${GPU_ARCH},cuda${CUDA_RUNTIME} -D__MPI -D__CUDA -D__GPU_MPI")
fi

# Manula version
#declare -a    flags=("" "-D__CUDA" "-D__MPI" "-D__MPI -D__CUDA" "-D__MPI -D__CUDA -D__GPU_MPI")
#declare -a    flags=("" "-D__MPI")

# END OF CONFIGURATIONS

function succ {
  local success=$1
  if [ $success -ne 0 ]; then
    echo "${red} Flags '$flag' FAILED ${reset}"
    exit $success
  else
    echo "${green} Flags '$flag' PASSED ${reset}"
  fi
}

# FUN
red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

for flag in "${flags[@]}"
do
  export DFLAGMATRIX="$flag"
  echo "Building with: $DFLAGMATRIX"
  cd ..
  make clean && make -f Makefile.test 1>/dev/null
  succ $?
  cd tests
  make clean && make 1>/dev/null
  succ $?
  for file in *.x
  do
    if [[ $flag = *"MPI"* ]]; then
      echo "Running $MPIEXEC ./$file"
      $MPIEXEC ./$file
      succ $?
    else
      echo "Running ./$file"
      ./$file
      succ $?
    fi
  done
done

