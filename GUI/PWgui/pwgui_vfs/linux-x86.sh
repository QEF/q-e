#!/bin/sh

ln -sf Makefile.linux Makefile
make
zip      pwgui-linux-x86.zip pwgui
tar zcvf pwgui-linux-x86.tgz pwgui
