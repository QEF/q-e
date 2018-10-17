#!/bin/bash



#modules
module load pgi/2018
module load cray-hdf5-parallel/1.10.2.0


#step in
cd q-e

espresso_version="6.3"
arch="v100"

#some stuff
export CUDA_ROOT=/project/projectdirs/mpccc/tkurth/cuda

#install directory
export SOFTWAREPATH=$(pwd)/install
export INSTALLPATH=${SOFTWAREPATH}/${espresso_version}/${arch}

#clean everything up so that no libraries with wrong arch are around
export FC=pgf90
export F90=pgf90
export MPIF90=mpif90
export FCFLAGS="-g -O3 -mp -Mpreprocess"
export F90FLAGS="${FCFLAGS}"
export MPIf90FLAGS="${F90FLAGS}"
export CC=pgcc
export CFLAGS="-g -O3 -mp"
export LDFLAGS="${F90FLAGS} -lgomp"

#hdf5
export HDF5_DIR=/home/tkurth/anaconda3

#configure
make veryclean
./configure --prefix=${INSTALLPATH} \
            --enable-openmp \
            --enable-parallel \
            --with-cuda=${CUDA_ROOT} \
            --with-cuda-cc=70 \
            --with-cuda-runtime=9.2 \
            --with-hdf5=${HDF5_DIR}/lib

#some makefile hacking
sed -i 's|^HDF5_LIB =|HDF5_LIB = -L${HDF5_DIR}/lib -lhdf5|g' make.inc

#clean up
make clean

#build crap
make pwall
#make want
#make gipaw
#make d3q
#make epw

#install
make install

#copy pseudo-dir:
cp -r pseudo ${INSTALLPATH}/


#step out
cd ..
