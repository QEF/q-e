
# Library of pseudopotential code

This directory contains a library of pseudopotential-related code,
extracted from the Quantum ESPRESSO distribution. This library is
stand-alone, except for some calls to LAPACK routines, and can be
independently compiled (a suitable `../make.inc` file is needed)

Currently, it includes basic definitions of the UPF (Unified Pseudopotential
File) format and I/O operations on them. UPF specifications are here:
http://www.quantum-espresso.org/pseudopotentials/unified-pseudopotential-format

In addition to the `Libupf.a` library, two executable utilities are produced:

- `upfconv.x`, converting pseudopotentials in other formats into UPF:
   see `upfconv.x -h` for more

- `virtual_v2.x`, courtesy Jingyang Wang (jw598@cornell.edu), generates
   an averaged pseudopotential suitable for Virtual Crystal Approximation

A python script `fixfile.py` is also present, to remove undesired `&`
characters from UPF files that hinder their parsing by xml tools.
