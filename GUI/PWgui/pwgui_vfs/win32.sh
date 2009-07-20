#!/bin/sh

ln -sf make.win32 make.sys
make

version=
if test -f make.versions; then
    version=`grep PWGUI_VERSION make.versions | awk '{print $NF}'`-
fi

zip      pwgui-${version}win32.zip pwgui.exe
tar zcvf pwgui-${version}win32.tgz pwgui.exe
