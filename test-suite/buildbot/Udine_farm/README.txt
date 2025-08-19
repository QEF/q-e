# To Install
# HDF 5.10.4

wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.4/src/hdf5-1.10.4.tar.gz

./configure --prefix=/home/buildbot/local/hdf5-104-gcc102 --enable-fortran

make -j 4
make install

###########################################################################################
# SP 21/01/25 - This can be easily copy-pasted to be in the same ENV as Buildbot to debug #
###########################################################################################
# 1) QE-farmer_gcc102_openmpi404
export BLAS_LIBS=/home/buildbot/local/blas380-gcc102/lib64/libblas.a
export LD_LIBRARY_PATH=/home/buildbot/local/lapack390-gcc102/lib:/home/buildbot/local/blas380-gcc102/lib64:/home/buildbot/local/openmpi-404-gcc102/lib/:/home/buildbot/local/gcc-10.2.0/lib64:/home/buildbot/local/gcc-10.2.0/lib:/home/buildbot/local/gcc-10.2.0/lib32
export PATH=/home/buildbot/local/openmpi-404-gcc102/bin:/home/buildbot/local/gcc-10.2.0/bin:/home/buildbot/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin

# 2) QE-farmer_gcc102_openmpi404-nimage
export BLAS_LIBS=/home/buildbot/local/blas380-gcc102/lib64/libblas.a
export LD_LIBRARY_PATH=/home/buildbot/local/lapack390-gcc102/lib:/home/buildbot/local/blas380-gcc102/lib64:/home/buildbot/local/openmpi-404-gcc102/lib/:/home/buildbot/local/gcc-10.2.0/lib64:/home/buildbot/local/gcc-10.2.0/lib:/home/buildbot/local/gcc-10.2.0/lib32
export PATH=/home/buildbot/local/openmpi-404-gcc102/bin:/home/buildbot/local/gcc-10.2.0/bin:/home/buildbot/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin

# QE-farmer_gcc102_openmpi404_hdf5
export LD_LIBRARY_PATH='/home/buildbot/local/hdf5-104-gcc102/lib:/home/buildbot/local/lapack390-gcc102/lib:/home/buildbot/local/blas380-gcc102/lib64:/home/buildbot/local/openmpi-404-gcc102/lib/:/home/buildbot/local/gcc-10.2.0/lib64:/home/buildbot/local/gcc-10.2.0/lib:/home/buildbot/local/gcc-10.2.0/lib32'
README.txt
export PATH='/home/buildbot/local/hdf5-104-gcc102/bin:/home/buildbot/local/openmpi-404-gcc102/bin:/home/buildbot/local/gcc-10.2.0/bin:/home/buildbot/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin'

# 3) QE-farmer_gcc102_openmpi404_openmp
export BLAS_LIBS=/home/buildbot/local/blas380-gcc102/lib64/libblas.a
export LD_LIBRARY_PATH=/home/buildbot/local/lapack390-gcc102/lib:/home/buildbot/local/blas380-gcc102/lib64:/home/buildbot/local/openmpi-404-gcc102/lib/:/home/buildbot/local/gcc-10.2.0/lib64:/home/buildbot/local/gcc-10.2.0/lib:/home/buildbot/local/gcc-10.2.0/lib32
export PATH=/home/buildbot/local/openmpi-404-gcc102/bin:/home/buildbot/local/gcc-10.2.0/bin:/home/buildbot/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin

# 4) QE-farmer_gcc102_serial
unset MKLROOT
export F77=gfortran
export F90=gfortran
export FC=gfortran
export CC=gcc
unset  CPP
unset  MPICC_CC
unset  MPICH_CC
unset  MPICH_CCC
unset  MPICH_CPP
unset  MPICH_CXX
unset  MPICH_F77
unset  MPICH_F90
unset  MPICXX_CXX
unset  MPIF90_F90
unset  OMPI_FC
unset  I_MPI_F90
unset  CXX

export BLAS_LIBS=/home/buildbot/local/blas380-gcc102/lib64/libblas.a
export LD_LIBRARY_PATH=/home/buildbot/local/lapack390-gcc102/lib:/home/buildbot/local/blas380-gcc102/lib64:/home/buildbot/local/gcc-10.2.0/lib64:/home/buildbot/local/gcc-10.2.0/lib:/home/buildbot/local/gcc-10.2.0/lib32
export PATH=/home/buildbot/local/gcc-10.2.0/bin:/usr/local/bin:/usr/bin:/bin:/snap/bin

# 5) QE-farmer_gcc750_openmpi316-dbg
export BLAS_LIBS=/home/buildbot/local/blas380-gcc102/lib64/libblas.a
export PATH=/home/buildbot/local/openmpi-316-gcc750/bin:/home/buildbot/local/gcc-750/bin:/home/buildbot/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin
export LD_LIBRARY_PATH=/home/buildbot/local/lapack390-gcc102/lib:/home/buildbot/local/blas380-gcc102/lib64:/home/buildbot/local/openmpi-316-gcc750/lib/:/home/buildbot/local/gcc-750/lib64:/home/buildbot/local/gcc-750/lib:/home/buildbot/local/gcc-750/lib32

# 6) QE-farmer_intel19_mvapich23
export LD_LIBRARY_PATH=/home/buildbot/local/mvapich2-23-intel19/lib:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/ipp/lib/intel64:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin
export PATH=/home/buildbot/local/mvapich2-23-intel19/bin:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/bin/intel64:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/bin:/home/buildbot/local/intel20/debugger_2020.4.304/gdb/intel64/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
export MKLROOT=/opt/intel/compilers_and_libraries_2020.4.304/linux/mkl
export MV2_SMP_USE_CMA=0

# 7) QE-farmer_intel20_impi
export MPICH_F90=ifort
export CC=icc
export MPICH_CC=icc
export FC=ifort
export F90=ifort
export F77=ifort
export MPICH_F77=ifort
export CPP='icc -E'
export MPICH_CPP='icc -E'
export CXX=icpc
export MPICH_CCC=icpc
export MPICH_CXX=icpc
export MPICC_CC=icc
export I_MPI_SHM_LMT=shm
export MPIF90_F90=mpiifort
export OMPI_FC=ifort
export I_MPI_F90=ifort
export MPICXX_CXX=icpc
export LD_LIBRARY_PATH=/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/compiler/lib/intel64_lin:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mpi/intel64/libfabric/lib:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mpi/intel64/lib/release:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mpi/intel64/lib:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/ipp/lib/intel64:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/tbb/lib/intel64/gcc4.8:/home/buildbot/local/intel20/debugger_2020/python/intel64/lib:/home/buildbot/local/intel20/debugger_2020/libipt/intel64/lib:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/daal/lib/intel64_lin:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/daal/../tbb/lib/intel64_lin/gcc4.4:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/daal/../tbb/lib/intel64_lin/gcc4.8
export PATH=/home/buildbot/local/intel20/intelpython3/bin:/home/buildbot/local/intel20/intelpython3/condabin:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/bin/intel64:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/bin:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mpi/intel64/libfabric/bin:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mpi/intel64/bin:/home/buildbot/local/intel20/debugger_2020/gdb/intel64/bin:/home/buildbot/buildbot/sandbox/bin:/home/buildbot/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/home/buildbot/local/intel20/parallel_studio_xe_2020.4.912/bin
export MKLROOT=/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mkl

# 8) QE-farmer_pgi2011_GPU
export MV2_SMP_USE_CMA=0
export CC=pgcc
export F90=nvfortran
export FC=nvfortran
export CUDA_HOME=/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/cuda/11.0
export OMPI_FC=nvfortran
export OPAL_PREFIX=/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/comm_libs/mpi
export MKLROOT=/opt/intel/compilers_and_libraries_2020.3.275/linux/mkl
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/cuda/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/math_libs/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/comm_libs/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/mpi/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/nccl/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/nvshmeme/lib:/opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/lib/intel64
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/cuda/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/compilers/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/20.11/comm_libs/mpi/bin:/home/buildbot/bin:/usr/local/bin:/usr/bin:/bin

# 9) QE-farmer_pgi2011_openmpi313
export MV2_SMP_USE_CMA=0
export CC=pgcc
export CXX=nvc++
export F77=nvfortran
export F90=nvfortran
export FC=nvfortran
export OMPI_FC=nvfortran
export MKLROOT=/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mkl
export LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/20.9/compilers/lib:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/20.9/compilers/lib:/home/buildbot/local/intel20/compilers_and_libraries_2020.4.304/linux/mkl/lib/intel64_lin
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/20.9/compilers/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/20.9/comm_libs/openmpi/openmpi-3.1.5/bin:/home/buildbot/.local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# 10) QE-farmer_pgi2107_GPU
export MV2_SMP_USE_CMA=0
export CC=pgcc
export F90=nvfortran
export FC=nvfortran
export CUDA_HOME=/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/cuda/11.4
export OMPI_FC=nvfortran
export OPAL_PREFIX=/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/comm_libs/mpi
export MKLROOT=/opt/intel/compilers_and_libraries_2020.3.275/linux/mkl
export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/cuda/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/compilers/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/math_libs/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/comm_libs/lib64:/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/compilers/mpi/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/compilers/nccl/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/compilers/nvshmeme/lib:/opt/intel/compilers_and_libraries_2020.3.275/linux/mkl/lib/intel64
export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/cuda/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/compilers/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/21.7/comm_libs/mpi/bin:/home/buildbot/bin:/usr/local/bin:/usr/bin:/bin























