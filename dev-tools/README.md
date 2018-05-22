This directory contains several tools that may be useful for developers

- `mem_counter`. A script that tracks all calls to `allocate` and `deallocate`,
   appending a call to subroutine `UtilXlib/mem_counter.f90`.
   Calls python script `mem_counter.py`, written by Pietro Bonfà (CINECA).
   `mem_counter -h` gives information on how to use it.
- `src-normal`. A script that "normalizes" the fortran syntax to QE style.
   Calls python script `src-normal.py`, written by Norbert Nemec.
   Usage: `src-normal file1.f90 [file2.f90 ...]` or `src-normal`
- `calltree.pl`
   A perl script, to be run from the root QE directory, producing in the
   standard output the tree of called routines 
- `callhtml.pl`
   As above, producing a html page with the tree of called routines 
- `release.sh`
   Script for packaging releases - obsolete, to be adapted to git
- utilities for PWgui:
  * `check_gui` (called via `Makefile`)
  * `diff_gui_help`
  * `guihelp.xsl`
  * `update_gui_help`
- utilities for helpdoc (see `README.helpdoc`):
  * `helpdoc`
  * `helpdoc.d`
  * `helpdoc.schema`
  * `input_xx.xsl`
- utilities for emacs_mode:
  * `gen-emacs-mode`
  * `gen-emacs-mode.tcl`
