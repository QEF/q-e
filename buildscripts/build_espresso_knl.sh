#!/bin/bash

#. ../environment.sh
export SOFTWAREPATH=$(pwd)/install
module load cray-hdf5/1.10.2.0
module load perftools-base
module load perftools-lite

#step in
cd q-e

espresso_version="6.3"
arch="knl"

#install directory
export INSTALLPATH=${SOFTWAREPATH}/${espresso_version}/${arch}

#clean everything up so that no libraries with wrong arch are around
compiler=intel
if [ "${compiler}" == "intel" ]; then
  export FC=ifort
  export F90=ifort
  export MPIF90=ftn
  export FCFLAGS="-g -O3 -qopenmp -mkl"
  export F90FLAGS="${FCFLAGS}"
  export MPIf90FLAGS="${F90FLAGS}"
  export CC=icc
  export CFLAGS="-g -O3 -qopenmp -mkl"
  export LDFLAGS="${F90FLAGS} -mkl"
fi

#configure
make veryclean
./configure --prefix=${INSTALLPATH} \
            --enable-openmp \
            --enable-parallel \
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
