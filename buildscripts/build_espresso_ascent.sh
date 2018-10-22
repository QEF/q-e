#!/bin/bash

#source environment
source env_ascent.sh

#install directory
export SOFTWAREPATH=$(pwd)/install
export INSTALLPATH=${SOFTWAREPATH}/${espresso_version}/${arch}

#step in
cd ../

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
rm -r ${INSTALLPATH}

#build crap
make -j 8 pwall
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
