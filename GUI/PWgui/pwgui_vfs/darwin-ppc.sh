#!/bin/sh

ln -sf make.darwin make.sys
make from_tarball
zip      pwgui-darwin-ppc.zip pwgui
tar zcvf pwgui-darwin-ppc.tgz pwgui
