#!/bin/sh

ln -sf Makefile.win32 Makefile
make
zip      pwgui-win32.zip pwgui.exe
tar zcvf pwgui-win32.tgz pwgui.exe
