#!/bin/sh

ln -sf make.darwin make.sys
make
zip      pwgui-darwin-ppc.zip pwgui
tar zcvf pwgui-darwin-ppc.tgz pwgui
