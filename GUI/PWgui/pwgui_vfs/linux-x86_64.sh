#!/bin/sh

ln -sf make.linux make.inc
make

# temporary workaround: 

if test ! -f pwgui; then
    # if transition from pwgui.kit to pwgui fails then we have a
    # non-optimal pwgui.kit, lets use it
    cp pwgui.kit pwgui
fi

version=
if test -f make.versions; then
    version=`grep PWGUI_VERSION make.versions | awk '{print $NF}'`-
fi

zip      pwgui-${version}linux-x86_64.zip pwgui
tar zcvf pwgui-${version}linux-x86_64.tgz pwgui
