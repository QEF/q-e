Quantum ESPRESSO GPU
====================

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

This repository contains the up to date GPU accelerated version of QuantumESPRESSO.

Installation
============

This version is tested against PGI compilers v. >= 17.4. The configure 
script checks the presence of a PGI compiler and of a few cuda libraries.
For this reason path pointing to cudatoolkit must be present in the
`LD_LIBRARY_PATH`.

A template for the configure command is:

```
./configure --with-cuda=XX --with-cuda-runtime=YY --with-cuda-cc=ZZ --enable-openmp [ --with-scalapack=no ]
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

From time to time PGI links to the wrong CUDA libraries anf fails reporting
a problem in `cusolver` missing `GOmp` (GNU Openmp). The solution to this
problem is removing cudatoolkit from the `LD_LIBRARY_PATH` before compiling.

Serial compilation is also supported.

Active branches
===============

There are currently two active branches:

* gpu_develop
* gpu_develop_multifft

Both tese branches are aligned with the develop branch of `QEF/q-e`.
The difference between the two branches is only in the FFTXlib library.
In the second case a 1D+2D algorithm is used for both the CPU and the GPU
subroutines. This makes the code considerably faster, especially on the
GPU side (generally a factor 2). This second branch should be considered
experimental.


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

The current GPU version passes all 186 tests with both parallel and serial
compilation. The testing suite should only be used to check the correctness of `pw.x`.
Therefore only `make run-tests-pw-parallel` and `make run-tests-pw-serial`
should be used.

Naming conventions
==================

Variables allocated on the device must end with `_d`.
Subroutines and functions replicating an algorithm on the GPU must end with `_gpu`.
Modules must end with `_gpum`.
Files with duplicated source code must end with `_gpu.f90`.

Porting functionalities
=======================

PW functionalities are ported to GPU by duplicating the subroutines and
the functions that operate on CPU variables.
The number of arguments should not change but input and output data may
be referring to device variables when applicable.

Bifurcations in code flow happen at runtime with commands similar to

```
use control_flags, only : use_gpu
[...]
if (use_gpu) then
   call subroutine_gpu(arg_d)
else
   call subroutine(arg)
end if
```

At each bifurcation point it should be possible to remove the call to the
accelerated routine without breaking the code. Note however that calling
both the CPU and the GPU version of a subroutine in the same place may
break the code execution.


Memory management
=================

[ DISCLAIMER STARTS ]
What described below is not the method that will be integrated
in the final release. Nonetheless it happens to be a good approach for:

1) simplify the alignment of this fork with the main repository,
2) debugging,
3) tracing evolution of memory paths as the CPU version evolves,
4) (in the future) report on a the set of global variables that should be 
   kept to guarantee a certain speedup.

For example, this simplified the integration of the changes that took
place to modernize the I/O.
[ DISCLAIMER ENDS ]


Global GPU data are tightly linked to global CPU data. One cannot allocate
global variables on the GPU manually. The global GPU variables follow the
allocation and deallocation of the CPU ones. This is an automatic mechanism
enforced by the managed memory system. In what follows, I will refer to 
duplicated GPU variables as "duplicated variable" and to the equivalent
CPU variable as "parent variable".

Global variables in modules are synchronized through calls to subroutines
named `using_xxx` and `using_xxx_d` with `xxx` being the name of the variable
in the module globally accessed by multiple subroutines.
This function accepts one argument that replicates the role of the `intent`
attribute.

Acceptable values are:
```
0: variable will only be read (equal to intent in)
1: variable will be read and written (equal to intent inout)
2: variable will be only (entirely) updated (equal to intent out).
```

Function and subroutine calls having global variables in their argument
should be guarded by calls to `using_xxx` with the appropriate argument.
Obviously calls with argument 0 and 1 must always be prepended.


The actual allocation of a duplicated variable happens when `using_xxx_d`
is called and the parent variable is allocated.
Deallocation happens when `using_xxx_d(2)` is called and the CPU variable
is not allocated.
Data synchronization (done with synchronous copies, i.e. overloaded cudamemcpy)
happens when either the CPU or the GPU memory is found to be flagged
"out of date" by a previous call to `using_xxx(1)` or `using_xxx(2)`
or `using_xxx_d(1)` or `using_xxx_d(2)`.

Calls to `using_xxx_d` should only happen in GPU function/subroutines.
This rule can be avoided if the call is protected by ifdefs.
This is useful if you are lazy and a global variable is updated only a few times.
An example of this being g vectors that are set in a few places (at
initialization, after a scaling of the Hamiltonian etc) and are used
everywhere in the code.

Finally, there are global variables that are only updated with subroutines
residing inside the same module. The allocation and the update of the
duplicated counterpart becomes trivial and is simply done at the same time
as the CPU variable. At the time of writing this constitute an exception
to the general rule but it is actually the result of the efforts done in
the last year to modularize the code and is probably the correct method
to deal with duplicated data in the code.
