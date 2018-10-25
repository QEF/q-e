#!/bin/bash

#modules
#very important to load gcc, otherwise cpp will cause issues
module load gcc/6.4.0
module load cuda/9.2.148
module load pgi/18.7
module load hdf5/1.10.3

#some parameters
espresso_version="6.3"
arch="p9_v100"

#some stuff
#export CUDA_ROOT=${OLCF_PGI_ROOT}/linuxpower/2018/cuda/9.2
export CUDA_ROOT=${OLCF_CUDA_ROOT}
#very important because the crappy pgi module sets it to the wrong value
export CPATH=

#clean everything up so that no libraries with wrong arch are around
export FC=pgf90
export F90=pgf90
export MPIF90=mpif90
export FCFLAGS="-g -O3 -mp -Mpreprocess -D__PROFILE -D__GPU_MPI"
export F90FLAGS="${FCFLAGS}"
export MPIf90FLAGS="${F90FLAGS}"
#export CPPFLAGS=-D__PROFILE
export CC=pgcc
export CFLAGS="-g -O3 -mp -Mpreprocess -D__PROFILE -D__GPU_MPI"
export LDFLAGS="${F90FLAGS}"
export NVCCFLAGS="-lineinfo"

#hdf5
export HDF5_DIR=${OLCF_HDF5_ROOT}
