#!/bin/sh -x

# Run this as "./dev-tools/release.sh"

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

#
VERSION=4.1
#
ESPRESSO_DIR=espresso-$VERSION
GUI=PWgui-$VERSION

# BEWARE: 
# in order to build the .html and .txt documentation in Doc, 
# "tcl", "tcllib", "xsltproc" are needed
# in order to build the GUI tarball from CVS sources, 
# "latex2html" and "convert" (from Image-Magick) are needed

if test -d $ESPRESSO_DIR; then /bin/rm -rf $ESPRESSO_DIR; fi

# produce updated ChangeLogs

make log
mv ChangeLog Doc/ChangeLog-$VERSION
mv ChangeLog.html Doc/ChangeLog-$VERSION.html

# produce documentation

cd doc-def/
make all
make clean
cd ../

cd Doc/
pdflatex user_guide.tex
pdflatex developer_man.tex
cd ../

# package using Makefile

make tar
make tar-gui PWGUI_VERSION=$VERSION

# unpackage in directory with version

mkdir $ESPRESSO_DIR
cd $ESPRESSO_DIR 
tar -xzf ../espresso.tar.gz
tar -xzf ../$GUI.tgz
/bin/rm ../$GUI.tgz ../espresso.tar.gz
cd ..
tar -cvzf espresso-$VERSION.tar.gz  $ESPRESSO_DIR >  espresso-$VERSION.lst
tar -cvzf $GUI.tar.gz $ESPRESSO_DIR/$GUI >  $GUI.lst
echo "espresso-$VERSION.tar.gz and  $GUI.tar.gz saved in directory:" `pwd`
echo "List of files in espresso-$VERSION.lst and $GUI.lst"

