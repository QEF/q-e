#!/bin/sh -x

VERSION=2.1.3
GUI_VERSION=2.1.3
GUI=PWgui-$GUI_VERSION

TMPDIR=./espresso-$VERSION

# check that $TMPDIR does not exist

if test -d $TMPDIR ; then (echo "$TMPDIR exists, stopping"; exit) ; fi
mkdir $TMPDIR

# make a copy of everything in $TMPDIR

tar -cf - bin/ config* README* Make* make* \
                      install-sh install/ moduledep.sh License upftools/ \
                      include/ Doc/ Modules/ clib/ flib/ \
                      PW/ PP/ PH/ Gamma/ PWNC/ PWCOND/ D3/ Raman/ pwtools/ \
		      CPV/ FPMD/ atomic/ atomic_doc/ examples/ pseudo/ $GUI \
         | (cd $TMPDIR ; tar -xf -) 

tar -czf $TMPDIR/cp-$VERSION.tar.gz \
                            $TMPDIR/bin/  $TMPDIR/config* $TMPDIR/README* \
                            $TMPDIR/Make* $TMPDIR/make*   $TMPDIR/install-sh \
                            $TMPDIR/install/ $TMPDIR/moduledep.sh \
                            $TMPDIR/License  $TMPDIR/upftools/     \
                            $TMPDIR/include/ $TMPDIR/Doc/ $TMPDIR/Modules/ \
                            $TMPDIR/clib/ $TMPDIR/flib/ \
                            $TMPDIR/CPV/ 

tar -czf $TMPDIR/fpmd-$VERSION.tar.gz \
                            $TMPDIR/bin/  $TMPDIR/config* $TMPDIR/README* \
                            $TMPDIR/Make* $TMPDIR/make*   $TMPDIR/install-sh \
                            $TMPDIR/install/ $TMPDIR/moduledep.sh \
                            $TMPDIR/License  $TMPDIR/upftools/     \
                            $TMPDIR/include/ $TMPDIR/Doc/ $TMPDIR/Modules/ \
                            $TMPDIR/clib/ $TMPDIR/flib/ \
                            $TMPDIR/FPMD/

tar -czf $TMPDIR/$GUI.tar.gz $TMPDIR/$GUI

tar -czf $TMPDIR/pw_src-$VERSION.tar.gz \
                            $TMPDIR/bin/  $TMPDIR/config* $TMPDIR/README* \
                            $TMPDIR/Make* $TMPDIR/make*   $TMPDIR/install-sh \
                            $TMPDIR/install/ $TMPDIR/moduledep.sh \
                            $TMPDIR/License  $TMPDIR/upftools/     \
                            $TMPDIR/include/ $TMPDIR/Doc/ $TMPDIR/Modules/ \
                            $TMPDIR/clib/ $TMPDIR/flib/ \
                            $TMPDIR/PW/   $TMPDIR/PP/ $TMPDIR/PH/ \
                            $TMPDIR/Gamma/ $TMPDIR/PWNC/ $TMPDIR/PWCOND/ \
                            $TMPDIR/D3/   $TMPDIR/Raman/ $TMPDIR/pwtools/

tar -czf $TMPDIR/examples-$VERSION.tar.gz $TMPDIR/examples/ $TMPDIR/pseudo/

tar -czf $TMPDIR/espresso-$VERSION.tar.gz \
                            $TMPDIR/bin/  $TMPDIR/config* $TMPDIR/README* \
                            $TMPDIR/Make* $TMPDIR/make*   $TMPDIR/install-sh \
                            $TMPDIR/install/ $TMPDIR/moduledep.sh \
                            $TMPDIR/License  $TMPDIR/upftools/     \
                            $TMPDIR/include/ $TMPDIR/Doc/ $TMPDIR/Modules/ \
                            $TMPDIR/clib/ $TMPDIR/flib/ \
                            $TMPDIR/PW/   $TMPDIR/PP/ $TMPDIR/PH/ \
                            $TMPDIR/Gamma/ $TMPDIR/PWNC/ $TMPDIR/PWCOND/ \
                            $TMPDIR/D3/   $TMPDIR/Raman/ $TMPDIR/pwtools/ \
                            $TMPDIR/CPV/ $TMPDIR/FPMD/ \
		            $TMPDIR/atomic/ $TMPDIR/atomic_doc/ \
                            $TMPDIR/examples/ $TMPDIR/pseudo/ $TMPDIR/$GUI

cd $TMPDIR

cp README     README-$VERSION
cp Doc/ChangeLog  ChangeLog-$VERSION
cp Doc/BUGS       BUGS-$VERSION
cp Doc/manual.tex manual-$VERSION.tex
cp Doc/manual.pdf manual-$VERSION.pdf

scp README-$VERSION ChangeLog-$VERSION BUGS-$VERSION \
    manual-$VERSION.tex manual-$VERSION.pdf \
    espresso-$VERSION.tar.gz examples-$VERSION.tar.gz $GUI.tar.gz  \
    pw_src-$VERSION.tar.gz cp-$VERSION.tar.gz fpmd-$VERSION.tar.gz \
    cibs.sns.it:public_html/pw

cd ../
/bin/rm -r $TMPDIR
