#!/bin/sh -x

version=6.2
revision=13897
user=giannozz
tempdir=$HOME/tempdir

# Make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

# work in $tempdir/qe-$version
cd $tempdir

# Get the svn copy via tag
# svn checkout http://qeforge.qe-forge.org/svn/q-e/tags/QE-$version/espresso qe-$version
# -OR- get the svn copy via revision checkout
svn checkout -r$revision svn+ssh://$user@qeforge.qe-forge.org/svnroot/q-e/trunk/espresso qe-$version
cd qe-$version

# Following operations require make.inc and svn files
touch make.inc

# Generate version.f90
cd Modules/
make version
/bin/rm version.f90.in
cd ..

# Generate documentation - NOTA BENE:
# in order to build the .html and .txt documentation in Doc,
# "tcl", "tcllib", "xsltproc" are needed
# in order to build the .pdf files in Doc, "pdflatex" is needed
# in order to build html files for user guide and developer manual,
# "latex2html" and "convert" (from Image-Magick) are needed

make doc VERSION=$version

# Generate PWGUI
make tar-gui PWGUI_VERSION=$version
mv PWgui-$version.tgz ../qe-$version-PWgui.tar.gz

# Generate QE-modes (requires tcllib, emacs, texlive-upquote)
make tar-qe-modes VERSION=$version
mv QE-modes-$version.tar.gz ../qe-$version-emacs_modes.tar.gz

# Move out svn directories, unneeded files, packages not to be packaged
/bin/rm .??* TODO archive/plumed-1.3-qe.tar.gz
/bin/rm -rf .svn/ QHA/ PlotPhon/ West/ GIPAW/ GUI/

# Update reference outputs on test-suite
cd test-suite/
for file in */benchmark.out*; do
    file2=`echo $file | sed 's/SVN/$version'`
    mv $file $file2
done
cd ..

# Package test-suite
cp License test-suite/
tar -czvf ../qe-$version-test-suite.tar.gz test-suite
/bin/rm -rf test-suite

# Package example
tar -czvf ../qe-$version-examples.tar.gz License */Examples */examples */*/examples
/bin/rm -rf */Examples */examples */*/examples

# Package sources (in directory qe-$version)
cd ..
tar -czvf qe-$version.tar.gz qe-$version
cd qe-$version

# Prepare documentation for upload
cd Doc/
/bin/rm *.xml *.txt *.html
cp ../*/Doc/*.html ../*/Doc/*.txt .

cp ../PW/Doc/user_guide.pdf ./pw_user_guide.pdf
cp ../CPV/Doc/user_guide.pdf ./cp_user_guide.pdf
cp ../PP/Doc/user_guide.pdf ./pp_user_guide.pdf
cp ../PHonon/Doc/user_guide.pdf ./ph_user_guide.pdf
cp ../NEB/Doc/user_guide.pdf ./neb_user_guide.pdf
cp ../atomic/Doc/pseudo-gen.pdf ./pseudo-gen.pdf

cp -R ../PW/Doc/user_guide ./pw_user_guide
cp -R ../CPV/Doc/user_guide ./cp_user_guide
cp -R ../PP/Doc/user_guide ./pp_user_guide
cp -R ../PHonon/Doc/user_guide ./ph_user_guide
cp -R ../NEB/Doc/user_guide ./neb_user_guide
cp -R ../atomic/Doc/pseudo-gen ./pseudo-gen

echo Now copy "Docs" to QE website
echo scp -R Doc user@site/wp-content/uploads/Doc-$version
echo connect to the website and create/update symbolic link to "Doc-$version"
