#!/bin/sh

ln -sf make.linux make.sys
make

version=
if test -f make.versions; then
    version=`grep PWGUI_VERSION make.versions | awk '{print $NF}'`-
fi

zip      pwgui-${version}linux-x86.zip pwgui
tar zcvf pwgui-${version}linux-x86.tgz pwgui
