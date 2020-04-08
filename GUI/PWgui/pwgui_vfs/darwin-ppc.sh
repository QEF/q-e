#!/bin/sh

ln -sf make.darwin make.inc
make

version=
if test -f make.versions; then
    version=`grep PWGUI_VERSION make.versions | awk '{print $NF}'`-
fi

zip      pwgui-${version}darwin-ppc.zip pwgui
tar zcvf pwgui-${version}darwin-ppc.tgz pwgui
