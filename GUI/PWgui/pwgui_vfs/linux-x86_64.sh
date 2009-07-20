#!/bin/sh

ln -sf make.linux make.sys
make

# temporary workaround: 

if test ! -f pwgui; then
    # transition from pwgui.kit --> pwgui fails at the moment, which
    # means we have non-optimal pwgui.kit, lets us it
    cp pwgui.kit pwgui
fi

version=
if test -f make.versions; then
    version=`grep PWGUI_VERSION make.versions | awk '{print $NF}'`-
fi

zip      pwgui-${version}linux-x86_64.zip pwgui
tar zcvf pwgui-${version}linux-x86_64.tgz pwgui
