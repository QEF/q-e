@echo off

rem ***
rem *** This is a template Batch file for MS-Windows (edit to suit your needs)
rem ***
rem *** EDIT: set the PWGUI variable to point to correct directory !!!

set initdir=%HOMEPATH%
set PWGUI=C:\Tone\PWgui-0.6.1

cd %initdir%
start wish %PWGUI%\pwgui.tcl