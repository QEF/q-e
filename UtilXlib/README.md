UtilXlib
========

This library implements various basic tasks such as timing, tracing,
optimized memory accesses and MPI communicators.

The following pre-processor directives can be used to enable/disable some features:

* `__MPI` : activates MPI support.
* `__TRACE` : activates verbose output for debugging purposes
* `__CUDA` : activates CUDA Fortran based interfaces.
* `__GPU_MPI` : use CUDA aware MPI calls instead of standard sync-send-update method (experimental).

Testing
=======

Partial unit testing is available in the `tests` sub-directory. See the 
README in that directory for further information.
