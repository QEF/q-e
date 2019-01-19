#!/bin/bash

#some parameters
espresso_version="6.3"
arch="knl"
compiler="pgi"

#install directory
export SOFTWAREPATH=$(pwd)/install
export INSTALLPATH=${SOFTWAREPATH}/${espresso_version}/${arch}/${compiler}

#step in
cd ../

#clean everything up so that no libraries with wrong arch are around
if [ "${compiler}" == "intel" ]; then
  module load intel
  fc=ifort
  f90=ifort
  mpif90=ftn
  fcflags="-g -O3 -qopenmp -mkl"
  f90flags="${fcflags}"
  mpif90flags="${f90flags}"
  cc=icc
  cflags=="-g -O3 -qopenmp -mkl"
  ldflags="${f90flags} -mkl"
fi
if [ "${compiler}" == "pgi" ]; then
  module load gcc
  module load pgi
  fc=pgf90
  f90=pgf90
  mpif90=pgf90
  fcflags="-O0 -mp -Mpreprocess -pgf90libs -I${MKLROOT}/include"
  f90flags="${fcflags}"
  mpif90flags="${f90flags}"
  cc=pgcc
  cflags="-O0 -mp -I${MKLROOT}/include"
  ldflags="${fcflags}"
fi

#configure
make veryclean
FC=${fc} F90=${f90} MPIF90=${mpif90} CC=${cc} LDFLAGS=${ldflags} ./configure --prefix=${INSTALLPATH} --enable-openmp --enable-parallel

#some makefile hacking
#sed -i 's|^HDF5_LIB =|HDF5_LIB = -L${HDF5_DIR}/lib -lhdf5|g' make.inc
if [ "${compiler}" == "intel" ]; then
    sed -r -i 's|^F90FLAGS\s{0,}=|F90FLAGS = -xMIC-AVX512 |g' make.inc
    sed -r -i 's|^FFLAGS\s{0,}=|FFLAGS = -xMIC-AVX512 |g' make.inc
    sed -r -i 's|^CFLAGS\s{0,}=|CFLAGS = -xMIC-AVX512 |g' make.inc
fi
if [ "${compiler}" == "pgi" ]; then
    sed -r -i 's|^F90FLAGS\s{0,}=|F90FLAGS = -tp=knl |g' make.inc
    sed -r -i 's|^FFLAGS\s{0,}=|FFLAGS = -tp=knl |g' make.inc
    sed -r -i 's|^CFLAGS\s{0,}=|CFLAGS = -tp=knl |g' make.inc
fi

#clean up
make clean
rm -rf ${INSTALLPATH}

#build crap
make pw
#make want
#make gipaw
#make d3q
#make epw

#install
mkdir -p ${INSTALLPATH}/bin
cp -r PW/src/*.x ${INSTALLPATH}/bin/

#copy pseudo-dir:
cp -r pseudo ${INSTALLPATH}/


#step out
cd ..
