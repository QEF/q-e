@echo off

rem ***
rem *** This is a template Batch file for MS-Windows (edit to suit your needs)
rem ***
rem *** EDIT: set the GUIB variable to point to the correct directory !!!

set GUIB=C:\Tone\Guib

start wish %GUIB%\guib.tcl %1
