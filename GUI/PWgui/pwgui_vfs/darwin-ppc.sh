#!/bin/sh

ln -sf Makefile.darwin Makefile
make
zip      pwgui-darwin-ppc.zip pwgui
tar zcvf pwgui-darwin-ppc.tgz pwgui
