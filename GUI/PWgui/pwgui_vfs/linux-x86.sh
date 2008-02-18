#!/bin/sh

ln -sf make.linux make.sys
make
zip      pwgui-linux-x86.zip pwgui
tar zcvf pwgui-linux-x86.tgz pwgui
