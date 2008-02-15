#!/bin/sh

ln -sf make.linux make.sys
make from_tarball
zip      pwgui-linux-x86_64.zip pwgui
tar zcvf pwgui-linux-x86_64.tgz pwgui
