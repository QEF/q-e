# Dev Tools

This directory contains several tools that may be useful for developers

- `mem_counter`. A script that tracks all calls to `allocate` and `deallocate`,
   appending a call to subroutine `UtilXlib/mem_counter.f90`.
   Calls python script `mem_counter.py`, written by Pietro Bonfà (CINECA)
   and improved by Samuel Poncé. 
   `mem_counter -h` gives information on how to use it.
   BEWARE: you may still need to manually edit some files that do not compile.
-  `mem_analyse.py` is a python script, by Samuel Poncé, locating memory leaks.
   See the script header for directions on how to use it.
- `src-normal`. A script that "normalizes" the fortran syntax to QE style
   (see below). Calls python script `src-normal.py`, written by Norbert Nemec.

   Usage: `src-normal file1.f90 [file2.f90 ...]` or `src-normal`
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

## Obsolescent utilities
- GPU utilities by Pietro Bonfà:
  * `get_device_props.py`
  * `device_props.c`
- Other utilities:
  * `calltree.pl`
   A perl script, to be run from the root QE directory, producing in the
   standard output the tree of called routines
  * `callhtml.pl`
   As above, producing a html page with the tree of called routines

## Coding style
These are some basic rules for Fortran codes enforced by `src_normal`:
* Use spaces for indentation instead of tabs (tab width 8 characters).
* Trailing whitespaces at the end the line should be removed.
* Normalize multiword keywords (e.g. END DO).
* Use capitalize version of the intrisic keywords (IF, DO, SUBROUTINE, etc.).
* Use the newest version of the comparison operators (==, >, etc.) instead of the old one (.eq., .gt., etc.)
