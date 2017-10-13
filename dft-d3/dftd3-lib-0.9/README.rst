=====
DFTD3
=====

This is a repackaged version of the `DFTD3 program
<http://www.thch.uni-bonn.de/tc/index.php?section=downloads&subsection=getd3>`_
by S. Grimme and his coworkers.

The original program (V3.1 Rev 1) was downloaded at 2016-04-03. It has been
converted to free format and encapsulated into modules. The source has been
split into two parts:

* A library with the core functionality. This can be directly used by third
  party applications wishing to calculate dispersion with the DFT-D3
  approach.
  
* Additional extensions which are necessary for the command line tool DFTD3 and
  the command line tool itself.

* Updated dftd3 code to include refitted/modified zero- and BJ-damped D3
  versions of Sherrill and coworkers (-bjm and -zerom)
  (Functionality corresponds to V3.2 Rev0)

Compilation
===========

Edit the file `make.arch` to reflect your compiler and linker. Then you can
issue one of the following commands:

* ``make lib``: to build the library `libdftd3.a` and the necessary
  module files (`*.mod`) in the directory `lib/`.

* ``make dftd3``: to build the executable `dftd3` in the directory `prg/`.

* ``make testapi``: to build a simple tester for the library (`testapi`) in the
  directory `test/`. The source code of this tester demonstrates how the library
  can be used by third party codes.

If you just issue ``make``, all three targets will be compiled.


Credits
=======

When using the library or the dftd3 tool, please cite:

  S. Grimme, J. Antony, S. Ehrlich and H. Krieg
  J. Chem. Phys, 132 (2010), 154104.
 
If BJ-damping is used 

  S. Grimme, S. Ehrlich and L. Goerigk
  J. Comput. Chem, 32 (2011), 1456-1465.

should be cited as well.


License
=======

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 1, or (at your option) any later version.
