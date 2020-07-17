
# Library of pseudopotential code

This directory contains a library of pseudopotential-related code,
extracted from the Quantum ESPRESSO distribution. This library is
almost stand-alone, except for the usage of FoX modules and routines,
for some calls to LAPACK routines, and for the need to include a
suitable `../make.inc` file in Makefile. Other than this, it can be
independently compiled.

Currently, it includes basic definitions of the UPF (Unified Pseudopotential
File) format and I/O operations on them. UPF specifications are here:
http://www.quantum-espresso.org/pseudopotentials/unified-pseudopotential-format

In addition to the `libupf.a` library, two executable utilities are produced:

- `upfconv.x`, converting pseudopotentials in other formats into UPF:
   see `upfconv.x -h` for more

- `virtual_v2.x`, courtesy Jingyang Wang (jw598@cornell.edu), generates
   an averaged pseudopotential suitable for Virtual Crystal Approximation

A python script `fixfile.py` is also present, to remove undesired `&`
characters from UPF files that hinder their parsing by xml tools.

