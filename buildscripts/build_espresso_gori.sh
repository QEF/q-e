#!/bin/bash

#source gori env
. env_gori.sh

#important variables
espresso_version="6.3"
arch="v100"
compiler="pgi"
#scalapackflags="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_pgi_thread -lmkl_core -lmkl_blacs_openmpi_lp64 -pgf90libs -mp -lpthread -lm -ldl"
export SOFTWAREPATH=$(pwd)/install

#step in
cd ../

#install directory
export INSTALLPATH=${SOFTWAREPATH}/${espresso_version}/${arch}/${compiler}

#clean everything up so that no libraries with wrong arch are around
fc=pgf90
fcflags="-g -O0 -mp -Mpreprocess -pgf90libs -I${MKLROOT}/include"
f90=pgf90
f90flags="${FCFLAGS}"
mpif90=mpif90
mpif90flags="${F90FLAGS}"
cc=pgcc
cflags="-g -O0 -mp -I${MKLROOT}/include"
ldflags="${F90FLAGS}"

#ompi stuff
export OMPI_CC=pgcc
export OMPI_FC=pgf90

#configure
make veryclean
FC=${fc} F90=${f90} MPIF90=${mpif90} CC=${cc} \
./configure --prefix=${INSTALLPATH} \
            --enable-openmp \
            --with-cuda=${CUDA_ROOT} \
            --with-cuda-cc=70 \
            --with-cuda-runtime=10.0 \
            --with-scalapack=no

#some makefile hacking
#sed -i 's|^HDF5_LIB =|HDF5_LIB = -L${HDF5_DIR}/lib -lhdf5|g' make.inc
#sed -i "s|^MANUAL_DFLAGS  =|MANUAL_DFLAGS  = -D__SCALAPACK|g" make.inc
#sed -i "s|^BLAS_LIBS      = ./*$|BLAS_LIBS      = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64  -lmkl_core -lmkl_pgi_thread|g" make.inc
sed -i "s|-D__DFTI|-D__FFTW|g" make.inc


#clean up
make clean

#build crap
make -j8 pw
#make want
#make gipaw
#make d3q
#make epw

#install
mkdir -p ${INSTALLPATH}/bin
cp PW/src/*.x ${INSTALLPATH}/bin/

#copy pseudo-dir:
cp -r pseudo ${INSTALLPATH}/

#step out
cd ..
