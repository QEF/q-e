#/bin/bash

#modules required
module load gcc/6.4.0
module load cuda/9.2.148
module load pgi/18.7
module load hdf5/1.10.3

#extract arguments
scratchdir=${1}

export OMP_NUM_THREADS=1

#extract some parameters
grank=$PMIX_RANK
lrank=$(($PMIX_RANK%6))

#step in
pushd ${scratchdir}

#nvprof runs on the first lrank (but OK to run on all lrank, for which each lrank would generate a .nvvp file)
if [[ ${lrank} == 0 ]];then
    APP="nvprof -o micro_%h_%p.nvvp ./pw.x -in micro.in -nbgrp ${nnodes}"  
else
    APP="./pw.x -in small.in -nbgrp ${nnodes}"
fi

use_hyperthreading=0  #0 not use ht; 1 use ht

if [[ $hyperthreading == 0 ]]
then
    case ${lrank} in
    [0])
    export PAMI_IBV_DEVICE_NAME=mlx5_0:1
    numactl --physcpubind=0,4,8,12,16,20,24 --membind=0 $APP
      ;;
    [1])
    export PAMI_IBV_DEVICE_NAME=mlx5_1:1
    numactl --physcpubind=28,32,36,40,44,48,52 --membind=0 $APP
      ;;
    [2])
    export PAMI_IBV_DEVICE_NAME=mlx5_0:1
    numactl --physcpubind=56,60,64,68,72,76,80 --membind=0 $APP
      ;;
    [3])
    export PAMI_IBV_DEVICE_NAME=mlx5_3:1
    numactl --physcpubind=88,92,96,100,104,108,112 --membind=8 $APP
      ;;
    [4])
    export PAMI_IBV_DEVICE_NAME=mlx5_2:1
    numactl --physcpubind=116,120,124,128,132,136,140 --membind=8 $APP
      ;;
    [5])
    export PAMI_IBV_DEVICE_NAME=mlx5_3:1
    numactl --physcpubind=144,148,152,156,160,164,168 --membind=8 $APP
      ;;
    esac
else 
    case ${lrank} in
    [0])
    export PAMI_IBV_DEVICE_NAME=mlx5_0:1
    numactl --physcpubind=0-27 --membind=0 $APP
      ;;
    [1])
    export PAMI_IBV_DEVICE_NAME=mlx5_1:1
    numactl --physcpubind=28-55 --membind=0 $APP
      ;;
    [2])
    export PAMI_IBV_DEVICE_NAME=mlx5_0:1
    numactl --physcpubind=56-83 --membind=0 $APP
      ;;
    [3])
    export PAMI_IBV_DEVICE_NAME=mlx5_3:1
    numactl --physcpubind=88-115 --membind=8 $APP
      ;;
    [4])
    export PAMI_IBV_DEVICE_NAME=mlx5_2:1
    numactl --physcpubind=116-143 --membind=8 $APP
      ;;
    [5])
    export PAMI_IBV_DEVICE_NAME=mlx5_3:1
    numactl --physcpubind=144-171 --membind=8 $APP
    ;;
    esac
fi
