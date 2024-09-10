Quantum ESPRESSO GPU
====================

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

This repository also contains the GPU-accelerated version of Quantum ESPRESSO.

Installation
============

This version requires the nvfortran (previously PGI) compiler from the
NVidia HPC SDK, v.21.7 or later (freely downloadable from NVidia). 
Earlier versions may or may not work and are no longer supported. 
You are advised to use the most recent version of NVidia software you can find. 

For compilation using CMake, see GitLab.com/QEF/q-e/-/wikis/Developers/CMake-build-system. For compilation using `configure`, see the User Guide in Doc/.
The `configure` script checks for the presence of the nvfortran compiler and 
of a few cuda libraries. For this reason the path pointing to the cuda toolkit
must be present in `LD_LIBRARY_PATH`. A template for the configure command is:

```
./configure --with-cuda=XX --with-cuda-runtime=YY --with-cuda-cc=ZZ --enable-openmp [ --with-scalapack=no ][ --with-cuda-mpi=yes ]
```

where `XX` is the location of the CUDA Toolkit (in HPC environments is 
typically `$NVHPC_CUDA_HOME` or `$CUDA_HOME`), `YY` is the version of 
the cuda toolkit and `ZZ` is the compute capability of the card. You can get 
those numbers from command `nvaccelinfo`, if you have a properly configured HPC SDK:
```
$ nvaccelinfo | grep -e 'Target' -e 'Driver'
CUDA Driver Version:           11000
Default Target:                cc70
...
```
The version is returned as (1000 major + 10 minor). For example, CUDA 11.0
is represented by 11000. For the above case, configure QE with:
```
./configure --with-cuda=$CUDA_HOME --with-cuda-cc=70 --with-cuda-runtime=11.0
```
One can also use command `nvidia-smi`: for two GPUs with cc70,
```
$ nvidia-smi --query-gpu=compute_cap --format=csv
7.0
7.0
```

Enabling faster communications between GPUs, via NVlink or Infiniband RDMA,
is essential for optimal performance. If your MPI library is built to be
CUDA-aware, then enable `--with-cuda-mpi=yes` (default: no). 

Serial (no MPI) compilation is also supported: use `--disable-parallel`.

Option --with-openacc is no longer honored: OpenACC is always needed.
It is generally a good idea to disable Scalapack when running small test
cases since the serial GPU eigensolver outperforms the parallel CPU
eigensolver in many circumstances.

From time to time PGI links to the wrong CUDA libraries and fails reporting a 
problem in `cusolver` missing `GOmp` (GNU Openmp). This problem can be solved
by removing the cuda toolkit from the `LD_LIBRARY_PATH` before compiling.

Execution
=========

By default, GPU support is active. The following message will appear at
the beginning of the output

```
     GPU acceleration is ACTIVE.
```

The current GPU version passes all tests with both parallel and serial 
compilation.
