Quantum ESPRESSO GPU
====================

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

This repository also contains the GPU-accelerated version of Quantum ESPRESSO.

Installation
============

This version requires the nvfortran (previously PGI) compiler from the
freely available NVidia HPC SDK. You are advised to use the most recent
version of NVidia software you can find. While any version later than 17.4
should work, many glitches are known to exist in older versions. 
The `configure` script checks for the presence of the nvfortran compiler and 
of a few cuda libraries. For this reason the path pointing to the cuda toolkit
must be present in `LD_LIBRARY_PATH`.

A template for the configure command is:

```
./configure --with-cuda=XX --with-cuda-runtime=YY --with-cuda-cc=ZZ --enable-openmp [--enable-openacc] [ --with-scalapack=no ]
```

where `XX` is the location of the CUDA Toolkit (in HPC environments is 
generally `$CUDA_HOME`), `YY` is the version of the cuda toolkit and `ZZ`
is the compute capability of the card. You can get those numbers from
command `nvaccelinfo`, if you have a properly configured HPC SDK:
```
$ nvaccelinfo | grep -e 'Target' -e 'Driver'
CUDA Driver Version:           11000
Default Target:                cc70
...
```
The version is returned as (1000 major + 10 minor). For example, CUDA 9.2 
would be represented by 9020. For the above case, configure QE with:
```
./configure --with-cuda=$CUDA_HOME --with-cuda-cc=70 --with-cuda-runtime=11.0
```
Alternatively, you may use the (deprecated) tool `get_device_props.py` in
directory `dev-tools/`.

It is generally a good idea to disable Scalapack when running small test
cases since the serial GPU eigensolver outperforms the parallel CPU
eigensolver in many circumstances.

From time to time PGI links to the wrong CUDA libraries and fails reporting a 
problem in `cusolver` missing `GOmp` (GNU Openmp). This problem can be solved
by removing the cuda toolkit from the `LD_LIBRARY_PATH` before compiling.

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
