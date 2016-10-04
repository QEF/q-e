#!/bin/sh -x

#tempdir=$HOME/Downloads
tempdir=/tmp
version=6.0
revision=13079

# make sure there is no locale setting creating unneeded differences.
#LC_ALL=C
#export LC_ALL

# get the svn copy via tag
svn checkout http://qeforge.qe-forge.org/svn/q-e/tags/QE-$version/espresso qe-$version

# -OR- get the svn copy via revision checkout
svn checkout -r$revision svn+ssh://<...>@qeforge.qe-forge.org/svnroot/q-e/trunk/espresso qe-$version

cd qe-$version

# *** manual edit Makefile  ***
# - Update PWgui
# - disable Doc distclean target

# *** manual edit install/plugins_makefile  ***
# - uncomment 'examples' target
# - uncomment 'uncompress-examples' target
# - uncomment 'examples_distclean' target

# Manual edit "userconfig.tmp" and "ENVIRONMENT"
# - change 'SVN' to $revision
# - change 'REFERENCE_VERSION' to $revision

# generate version.f90 (requires svn files)
# save version.f90 (make veryclean removes it)
touch make.inc
cd Modules
make version.f90
mv version.f90 ..
cd ..

# remove all .svn directories, clean
find . -type d -name .svn -exec /bin/rm -rf {} \;
make veryclean
rm archive/plumed-1.3-qe.tar.gz archive/PLUMED-latest.tar.gz

# restore version.f90
mv version.f90 Modules/
cp Modules/version.f90 Modules/version.f90.in
chmod -x install/update_version

# generate documentation - NOTA BENE:
# in order to build the .html and .txt documentation in Doc,
# "tcl", "tcllib", "xsltproc" are needed
# in order to build the .pdf files in Doc, "pdflatex" is needed
# in order to build html files for user guide and developer manual,
# "latex2html" and "convert" (from Image-Magick) are needed
touch make.inc
make doc VERSION=$version

# generate PWGUI
make tar-gui PWGUI_VERSION=$version
tar -xzvf PWgui-$version.tgz
/bin/rm PWgui-$version.tgz

# generate QE-modes (requires tcllib and emacs)
make tar-qe-modes VERSION=$version
mv QE-modes-$version.tar.gz ../qe-$version-emacs_modes.tar.gz

# Updating reference outputs on test-suite
cd test-suite
find . -name benchmark.out* > list-SVN.txt
sed 's/SVN/'$version'/g' list-SVN.txt | grep -v svn  > list-$version.txt
paste -d " " list-SVN.txt list-$version.txt > ./STUFF-TO-RENAME.txt
IFS=$'\n'
for x in `cat ./STUFF-TO-RENAME.txt `
do
file_src=`echo $x | awk '{ print $1}'`
file_dst=`echo $x | awk '{ print $2}'`
mv ${file_src} ${file_dst}
done
rm ./STUFF-TO-RENAME.txt ./list-SVN.txt ./list-$version.txt
cp License test-suite/

cd ..

make distclean

# packacking test-suite
mv test-suite test-suite
tar -czvf ../qe-$version-test-suite.tar.gz test-suite

# Grouping Examples in the same directory and packacking them
mkdir Examples
cd  Examples
mkdir CPV PHonon NEB COUPLE PP PW XSpectra GWW EPW atomic PWCOND TDDFPT PWgui-$version
mkdir PP/simple_transport
mv ../TDDFPT/Examples/* TDDFPT/
mv ../PWgui-$version/examples/* PWgui-$version/
mv ../atomic/examples/* atomic/
mv ../EPW/examples/* EPW/
mv ../GWW/examples/* GWW/
mv ../XSpectra/examples/* XSpectra/
mv ../PW/examples/* PW/
mv ../PP/examples/* PP/
mv ../PP/simple_transport/examples/* PP/simple_transport/
mv ../COUPLE/examples/* COUPLE/
mv ../NEB/examples/* NEB/
mv ../CPV/examples/* CPV/
mv ../PWCOND/examples/* PWCOND/
mv ../PHonon/examples/* PHonon/
rm -rf ../TDDFPT/Examples ../CPV/examples ../PHonon/examples ../NEB/examples ../COUPLE/examples ../PP/examples ../PP/simple_transport/examples ../PW/examples ../PWgui-6.0/examples ../XSpectra/examples ../GWW/examples ../EPW/examples ../atomic/examples ../PWCOND/examples
cd ..
cp License Examples/

# Grouping Examples in the same directory and packacking them
tar -czvf ../qe-$version-examples.tar.gz Examples

cd ../

# core espresso
tar -czvf qe-$version.tar.gz \
qe-$version/COUPLE \
qe-$version/CPV \
qe-$version/Doc \
qe-$version/EPW \
qe-$version/FFTXlib \
qe-$version/GWW \
qe-$version/LAXlib \
qe-$version/LR_Modules \
qe-$version/License \
qe-$version/Makefile \
qe-$version/Modules \
qe-$version/NEB \
qe-$version/PHonon \
qe-$version/PP \
qe-$version/PW \
qe-$version/PWCOND \
qe-$version/README \
qe-$version/TDDFPT \
qe-$version/XSpectra \
qe-$version/archive \
qe-$version/atomic \
qe-$version/clib \
qe-$version/configure \
qe-$version/dev-tools \
qe-$version/environment_variables \
qe-$version/include \
qe-$version/install \
qe-$version/pseudo \
qe-$version/upftools \
qe-$version/PWgui-$version

cd qe-$version

# Preparing documentation for upload
cd Doc
rm *.xml *.txt *.html

for x in CPV PHonon NEB PP PW PWCOND TDDFPT atomic; \
  do cp ../$x/Doc/*.xml .; cp ../$x/Doc/*.html .; cp ../$x/Doc/*.txt .; done

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

# Copy "Docs" to QE website
scp -R Doc <...>@<...>/wp-content/uploads/Doc-$version

# Connect to the website and create/update symbolic link to "Doc-$version"
