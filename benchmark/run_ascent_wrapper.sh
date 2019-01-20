#/bin/bash

#modules required
module load gcc/6.4.0
module load cuda/9.2.148
module load pgi/18.7
module load hdf5/1.10.3

#hack some OMPI params
export OMPI_MCA_osc_pami_allow_thread_multiple=0
export PAMI_ENABLE_STRIPING=0
export PAMI_IBV_ENABLE_OOO_AR=0
export PAMI_IBV_QP_SERVICE_LEVEL=0

#enable QE GPU suppoert
export USEGPU=1
#export PGI_TERM='trace'

#OpenMP tweaks
export OMP_NUM_THREADS=1

#extract arguments
scratchdir=${1}

#extract somne parameters
grank=$PMIX_RANK
lrank=$(($PMIX_RANK%6))

#step in
pushd ${scratchdir}

NVPROF="nvprof -f -o micro_multifft.nvvp"
APP="./pw.x -in large.in -nbgrp 1"


case ${lrank} in
[0])
export CUDA_VISIBLE_DEVICES=0
export PAMI_IBV_DEVICE_NAME=mlx5_0:1
numactl -N 0 $APP
#numactl --physcpubind=0,4,8,12,16,20,24 --membind=0 $APP
#numactl --physcpubind=0-27 --membind=0 $APP
  ;;
[1])
export CUDA_VISIBLE_DEVICES=1
export PAMI_IBV_DEVICE_NAME=mlx5_1:1
numactl -N 0 $APP
#numactl --physcpubind=28,32,36,40,44,48,52 --membind=0 $APP
#numactl --physcpubind=28-55 --membind=0 $APP
  ;;
[2])
export CUDA_VISIBLE_DEVICES=3
export PAMI_IBV_DEVICE_NAME=mlx5_3:1
numactl -N 8 $APP
#numactl --physcpubind=88,92,96,100,104,108,112 --membind=8 $APP
#numactl --physcpubind=88-115 --membind=8 $APP
  ;;
[3])
export CUDA_VISIBLE_DEVICES=4
export PAMI_IBV_DEVICE_NAME=mlx5_2:1
numactl -N 8 $APP
#numactl --physcpubind=116,120,124,128,132,136,140 --membind=8 $APP
#numactl --physcpubind=116-143 --membind=8 $APP   
  ;;
[4])
export CUDA_VISIBLE_DEVICES=2
export PAMI_IBV_DEVICE_NAME=mlx5_0:1
numactl -N 0 $APP
#numactl --physcpubind=56,60,64,68,72,76,80 --membind=0 $APP
#numactl --physcpubind=56-83 --membind=0 $APP      
  ;;
[5])
export CUDA_VISIBLE_DEVICES=5
export PAMI_IBV_DEVICE_NAME=mlx5_3:1
numactl -N 8 $APP
#numactl --physcpubind=144,148,152,156,160,164,168 --membind=8 $APP
#numactl --physcpubind=144-171 --membind=8 $APP
  ;;
esac
