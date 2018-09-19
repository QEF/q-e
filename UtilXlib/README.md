UtilXlib
========

This library implements various basic tasks such as timing, tracing,
optimized memory accesses and MPI communicators.

The following pre-processor directives can be used to enable/disable some
features:

* `__MPI` : activates MPI support.
* `__TRACE` : activates verbose output for debugging purposes
* `__CUDA` : activates CUDA Fortran based interfaces.
* `__GPU_MPI` : use CUDA aware MPI calls instead of standard sync-send-update method (experimental).

CUDA specific notes
===================

All calls to message passing interfaces are synchronous with respect to
both MPI and CUDA streams. The code will synchronize the device before
starting the communication, also in those cases where communication
may be avoided (for example in serial version).
A different behaviour may be observed when the default stream 
synchronization behaviour is overridden by the user (see `cudaStreamCreateWithFlags`).

Be careful when using CUDA-aware MPI. Some implementations are not
complete. The library will not check for the CUDA-aware MPI APIs during
the initialization, but may report failure codes during the execution.
If you encounter problems when adding the flag `__GPU_MPI` it might
be that the MPI library does not support some CUDA-aware APIs.


Known Issues
============
Owing to the use of the `source` option in data allocations,
PGI versions older than 17.10 may fail with arrays having initial index
different from 1.

Testing
=======

Partial unit testing is available in the `tests` sub-directory. See the 
README in that directory for further information.
