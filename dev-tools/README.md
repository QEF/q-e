# Dev Tools

This directory contains several tools that may be useful for developers

- `mem_counter`. A script that tracks all calls to `allocate` and `deallocate`,
   appending a call to subroutine `UtilXlib/mem_counter.f90`.
   Calls python script `mem_counter.py`, written by Pietro Bonfà (CINECA).
   `mem_counter -h` gives information on how to use it.
- `src-normal`. A script that "normalizes" the fortran syntax to [QE style](#style).
   Calls python script `src-normal.py`, written by Norbert Nemec.

   Usage: `src-normal file1.f90 [file2.f90 ...]` or `src-normal`

- `calltree.pl`
   A perl script, to be run from the root QE directory, producing in the
   standard output the tree of called routines
- `callhtml.pl`
   As above, producing a html page with the tree of called routines
- `release.sh`
   Script for packaging releases - obsolete, to be adapted to git
- Utilities for PWgui:
  * `check_gui` (called via `Makefile`)
  * `diff_gui_help`
  * `guihelp.xsl`
  * `update_gui_help`
- Utilities for helpdoc (see `README.helpdoc`):
  * `helpdoc`
  * `helpdoc.d`
  * `helpdoc.schema`
  * `input_xx.xsl`
- Utilities for emacs_mode:
  * `gen-emacs-mode`
  * `gen-emacs-mode.tcl`

## Coding style
This are some basic rules to follow when writting Fortran code.
* Use spaces for indentation instead of tabs (tab width 8 characters).
* Trailing whitespaces at the end the line should be removed.
* Normalize multiword keywords.
* Use capitalize version of the intrisic keywords (IF, DO, SUBROUTINE, etc.).
* Prefer the newest version of the comparison operators (==, >, etc.) instead of the old one (.eq., .gt., etc.)
