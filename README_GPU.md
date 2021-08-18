Quantum ESPRESSO GPU
====================

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

This repository also contains the GPU-accelerated version of Quantum ESPRESSO.

Installation
============

This version requires the nvfortran (previously PGI) compiler from the
freely available NVidia HPC SDK. You are adviced to use a recent version
of NVidia software. Any version later than 17.4 should work, but many glitches
are know to exist in older versions. 
The configure script checks for the presence of the nvfortran compiler and of 
a few cuda libraries.For this reason the path pointing to cudatoolkit must be
present in `LD_LIBRARY_PATH`.

A template for the configure command is:

```
./configure --with-cuda=XX --with-cuda-runtime=YY --with-cuda-cc=ZZ --enable-openmp [--enable-openacc] [ --with-scalapack=no ]
```

where `XX` is the location of the CUDA Toolkit (in HPC environments is 
generally `$CUDA_HOME`), `YY` is the version of the cuda toolkit and `ZZ`
is the compute capability of the card. 
If you have no idea what these numbers are you may give a try to the
automatic tool `get_device_props.py`. An example using Slurm is:

```
$ module load cuda
$ cd dev-tools
$ salloc -n1 -t1
[...]
salloc: Granted job allocation xxxx
$ srun python get_device_props.py
[...]
Compute capabilities for dev 0: 6.0
Compute capabilities for dev 1: 6.0
Compute capabilities for dev 2: 6.0
Compute capabilities for dev 3: 6.0

 If all compute capabilities match, configure QE with:
./configure --with-cuda=$CUDA_HOME --with-cuda-cc=60 --with-cuda-runtime=9.2
```

It is generally a good idea to disable Scalapack when running small test
cases since the serial GPU eigensolver can outperform the parallel CPU
eigensolver in many circumstances.

From time to time PGI links to the wrong CUDA libraries and fails reporting
a problem in `cusolver` missing `GOmp` (GNU Openmp). The solution to this
problem is removing cudatoolkit from the `LD_LIBRARY_PATH` before compiling.

Serial compilation is also supported.

Execution
=========

By default, GPU support is active. The following message will appear at
the beginning of the output

```
     GPU acceleration is ACTIVE.
```

GPU acceleration can be switched off by setting the following environment
variable:

```
$ export USEGPU=no
```


Testing
=======

The current GPU version passes all tests with both parallel and serial 
compilation.
