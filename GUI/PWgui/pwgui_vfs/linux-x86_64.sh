#!/bin/sh

ln -sf make.linux make.sys
make

# temporary workaround: 

if test ! -f pwgui; then
    # transition from pwgui.kit --> pwgui fails at the moment, which
    # means wem have non-optimal pwgui.kit, lets us it
    cp pwgui.kit pwgui
fi


zip      pwgui-linux-x86_64.zip pwgui
tar zcvf pwgui-linux-x86_64.tgz pwgui
