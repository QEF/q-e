#!/bin/sh -x

# BEWARE: 
# in order to build the GUI tarball from CVS sources the following
# software is needed:
#   1. pdflatex
#   2. convert (from Image-Magick)
#   3. latex2html

dir=pw-2-1
cd $HOME
if test -d $dir.save; then /bin/rm -r $dir.save; fi
if test -d $dir; then mv $dir $dir.save; fi
mkdir $dir

if test -d espresso ; then
    cd espresso
else
    echo "
   Oops. espresso/ does not exist. Aborting ...
"
    exit 1
fi

GUI_VERSION=`cat GUI/PWgui/VERSION`
GUI=PWgui-$GUI_VERSION
make veryclean
find . -type f -name *~ -exec /bin/rm {} \;
find . -type f -name .#* -exec /bin/rm {} \;
if test -f pw.tar.gz ; then /bin/rm pw.tar.gz ; fi
make tar tar-gui
cd  ../$dir

mkdir bin
tar -xzf ../O-sesame/pw.tar.gz
tar -xzf ../O-sesame/$GUI.tgz
find $GUI -name CVS -exec /bin/rm -r {} \;

tar -czf ../cp.tar.gz bin/ config* README* Make* make*           \
                      install-sh install/ moduledep.sh License upftools/ \
                      include/ Doc/ Modules/ clib/ flib/ \
                      CPV/

tar -czf ../fpmd.tar.gz bin/ config* README* Make* make*        \
                      install-sh install/ moduledep.sh License upftools/ \
                      include/ Doc/ Modules/ clib/ flib/ \
                      FPMD/

tar -czf ../$GUI.tar.gz $GUI

tar -czf ../pw.tar.gz bin/ config* README* Make* make* \
                      install-sh install/ moduledep.sh License upftools/ \
                      include/ Doc/ Modules/ clib/ flib/ \
                      PW/ PP/ PH/ Gamma/ PWNC/ PWCOND/ D3/ Raman/ pwtools/

tar -czf ../examples.tar.gz examples/ pseudo/

tar -czf ../espresso.tar.gz  bin/ config* README* Make* make* \
                      install-sh install/ moduledep.sh License upftools/ \
                      include/ Doc/ Modules/ clib/ flib/ \
                      PW/ PP/ PH/ Gamma/ PWNC/ PWCOND/ D3/ Raman/ pwtools/ \
		      CPV/ FPMD/ atomic/ atomic_doc/ examples/ pseudo/ $GUI

scp Doc/README Doc/ChangeLog Doc/BUGS Doc/manual.tex \
    Doc/*.png  Doc/manual.pdf ../pw.tar.gz ../espresso.tar.gz \
    ../examples.tar.gz ../$GUI.tar.gz ../fpmd.tar.gz ../cp.tar.gz \
    cibs:public_html/pw

