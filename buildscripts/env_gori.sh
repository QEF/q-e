#!/bin/bash

#unload crap
module purge
module load escori
#module unload altd
#module unload craype-haswell
#module unload cray-mpich
#load required modules
module load cuda/10.0
module load pgi/2018

#important paths and variables
mkldir=/global/common/cori/software/intel/compilers_and_libraries_2018.5.274/linux/mkl
mpipath=/global/homes/t/tkurth/src/openmpi_ucx/install_pgi/ompi

#source MKL
. ${mkldir}/bin/mklvars.sh intel64


#mpi stuff
#path hacking
export LD_LIBRARY_PATH=/opt/esslurm/lib64:${mpipath}/lib/:${LD_LIBRARY_PATH}
export PATH=${mpipath}/bin/:${PATH}

#openmpi env stuff
export OMPI_MCA_btl_openib_if_include=mlx5_0,mlx5_2,mlx5_4,mlx5_6
