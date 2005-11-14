#!/bin/sh -x

VERSION=3.0
TMPDIR=espresso-$VERSION

TARGET_MACHINE=cibs.sns.it:public_html/public/

GUI_VERSION=`cat GUI/PWgui/VERSION`
GUI=PWgui-$GUI_VERSION

# BEWARE: 
# in order to build the GUI tarball from CVS sources the following
# software is needed:
#   1. pdflatex
#   2. convert (from Image-Magick)
#   3. latex2html

if test -d $TMPDIR.save; then /bin/rm -r $TMPDIR.save; fi
if test -d $TMPDIR; then mv $TMPDIR $TMPDIR.save; fi
mkdir $TMPDIR
mkdir $TMPDIR/bin

make veryclean
find . -type f -name *~ -exec /bin/rm {} \;
find . -type f -name .#* -exec /bin/rm {} \;
if test -f espresso.tar.gz ; then /bin/rm espresso.tar.gz ; fi
##### Using Makefile:
make tar tar-gui
##### Alternatively:
#tar -czf espresso.tar.gz config* README* Make* make* \
#                         install-sh install/ moduledep.sh includedeps.sh License upftools/ \
#                         include/ Doc/ Modules/ clib/ flib/ \
#                         PW/ PP/ PH/ Gamma/ PWCOND/ D3/ pwtools/ \
#                         CPV/ atomic/ atomic_doc/ examples/ pseudo/
#tar -cf $GUI.tgz $GUI
##### End
cd  $TMPDIR
tar -xzf ../espresso.tar.gz
tar -xzf ../$GUI.tgz
find . -name CVS -type d -exec /bin/rm -r {} \;

cd ../

tar -czf $TMPDIR/cp-$VERSION.tar.gz \
                            $TMPDIR/bin/     $TMPDIR/config* $TMPDIR/README* \
                            $TMPDIR/Make*    $TMPDIR/make*   $TMPDIR/install-sh \
                            $TMPDIR/install/ $TMPDIR/moduledep.sh $TMPDIR/includedep.sh \
                            $TMPDIR/License  $TMPDIR/upftools/     \
                            $TMPDIR/include/ $TMPDIR/Doc/    $TMPDIR/Modules/ \
                            $TMPDIR/iotk/    $TMPDIR/clib/    $TMPDIR/flib/ \
                            $TMPDIR/CPV/ 

tar -czf $TMPDIR/$GUI.tar.gz $TMPDIR/$GUI

tar -czf $TMPDIR/pw_src-$VERSION.tar.gz \
                            $TMPDIR/bin/     $TMPDIR/config* $TMPDIR/README* \
                            $TMPDIR/Make*    $TMPDIR/make*   $TMPDIR/install-sh \
                            $TMPDIR/install/ $TMPDIR/moduledep.sh  $TMPDIR/includedep.sh \
                            $TMPDIR/License  $TMPDIR/upftools/     \
                            $TMPDIR/include/ $TMPDIR/Doc/ $TMPDIR/Modules/ \
                            $TMPDIR/iotk/    $TMPDIR/clib/ $TMPDIR/flib/ \
                            $TMPDIR/PW/      $TMPDIR/PP/ $TMPDIR/PH/ \
                            $TMPDIR/Gamma/   $TMPDIR/PWCOND/ \
                            $TMPDIR/D3/      $TMPDIR/pwtools/

tar -czf $TMPDIR/examples-$VERSION.tar.gz $TMPDIR/examples/ $TMPDIR/pseudo/


tar -czf $TMPDIR/espresso-$VERSION.tar.gz \
                            $TMPDIR/bin/     $TMPDIR/config* $TMPDIR/README* \
                            $TMPDIR/Make*    $TMPDIR/make*   $TMPDIR/install-sh \
                            $TMPDIR/install/ $TMPDIR/moduledep.sh $TMPDIR/includedep.sh \
                            $TMPDIR/License  $TMPDIR/upftools/     \
                            $TMPDIR/include/ $TMPDIR/Doc/    $TMPDIR/Modules/ \
                            $TMPDIR/iotk/    $TMPDIR/clib/   $TMPDIR/flib/ \
                            $TMPDIR/PW/      $TMPDIR/PP/ $TMPDIR/PH/ \
                            $TMPDIR/Gamma/   $TMPDIR/PWCOND/ \
                            $TMPDIR/D3/      $TMPDIR/pwtools/ \
                            $TMPDIR/CPV/  \
		            $TMPDIR/atomic/  $TMPDIR/atomic_doc/ \
                            $TMPDIR/pseudo/  $TMPDIR/examples/ $TMPDIR/$GUI
cd $TMPDIR

cp README         README-$VERSION
cp Doc/ChangeLog  ChangeLog-$VERSION
cp Doc/BUGS       BUGS-$VERSION
cp Doc/manual.tex manual-$VERSION.tex
cp Doc/manual.pdf manual-$VERSION.pdf

scp README-$VERSION ChangeLog-$VERSION BUGS-$VERSION \
    manual-$VERSION.tex manual-$VERSION.pdf \
    espresso-$VERSION.tar.gz examples-$VERSION.tar.gz $GUI.tar.gz  \
    pw_src-$VERSION.tar.gz cp-$VERSION.tar.gz \
    $TARGET_MACHINE

#cd ../
#/bin/rm -r $TMPDIR
