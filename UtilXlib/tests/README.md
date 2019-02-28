# UtilXlib testing suite

A set of tests for UtilXlib also showing the functionalities of the library.

In order to run the tests first run `./configure` in QE topdir and generate a
valid make.inc
After that edit the variables at the top of `compile_and_run_tests.sh` and run with:

    bash compile_and_run_tests.sh -sm[cn]

Options meaning:

 * -s: serial compilation and execution
 * -m: MPI compilation and execution
 * -c: CUDA-Fortran interface with data transfer done with memory on the *host*
 * -n: CUDA-Fortran interface with data transfer done with memory on the *device*

The tester module is part of the fortran_tester project avaiulable at:
https://github.com/pdebuyl/fortran_tester/archive/master.zip

TODO:

 * Get rid of seeds printed on screen. It's pointless.
