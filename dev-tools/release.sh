#!/bin/sh -x

#tempdir=$HOME/Downloads
tempdir=/tmp
version=6.0

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

mkdir $tempdir
cd $tempdir
/bin/rm -rf espresso/ qe-$version
# get the svn copy
svn checkout http://qeforge.qe-forge.org/svn/q-e/tags/QE-5.3.0/espresso
mv espresso/ qe-$version/

cd qe-$version

# generate version.f90 (requires svn files)
touch make.unc
cd Modules
make version.f90
# save version.f90 (make veryclean removes it)
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

# generate QE-modes

make tar-qe-modes VERSION=$version; # this creates a ready to use QE-modes-$version.tar.gz
# move the package one directory down from espresso-$version/
mv QE-modes-$version.tar.gz ..

cd ..

# Updating reference outputs on test-suite

find . -name benchmark.out* > list-SVN.txt
sed 's/SVN/$version/g' list-SVN.txt | grep -v svn  > list-$version.txt
paste -d " " list-SVN.txt list-$version.txt > ./STUFF-TO-RENAME.txt
IFS=$'\n'
for x in `cat ./STUFF-TO-RENAME.txt `
do
file_src=`echo $x | awk '{ print $1}'`
file_dst=`echo $x | awk '{ print $2}'`
mv ${file_src} ${file_dst}
done
rm ./STUFF-TO-RENAME.txt ./list-SVN.txt ./list-$version.txt

# core espresso

tar -czvf qe-$version.tar.gz qe-$version/archive \
                                   qe-$version/clib \
                                   qe-$version/configure \
                                   qe-$version/COUPLE \
                                   qe-$version/CPV \
                                   qe-$version/dev-tools \
                                   qe-$version/Doc \
                                   qe-$version/environment_variables \
                                   qe-$version/LAXlib \
                                   qe-$version/FFTXlib \
                                   qe-$version/Makefile \
                                   qe-$version/include \
                                   qe-$version/install \
                                   qe-$version/License \
                                   qe-$version/Modules \
                                   qe-$version/PP \
                                   qe-$version/pseudo \
                                   qe-$version/PW \
                                   qe-$version/README \
                                   qe-$version/upftools \
                                   qe-$version/NEB \
                                   qe-$version/PHonon \
                                   qe-$version/XSpectra \
                                   qe-$version/PWCOND \
                                   qe-$version/TDDFPT \
                                   qe-$version/atomic \
                                   qe-$version/gwl \
                                   qe-$version/atomic

# Packages, ready for automatic unpacking

cd qe-$version
tar -cvzf ../PWgui-$version.tar.gz    PWgui-$version
tar -czvf ../EPW-$version.tar.gz      EPW
tar -czvf ../GWW-$version.tar.gz      GWW
#tar -czvf ../GIPAW-$version.tar.gz    GIPAW
tar -czvf ../test-suite-$version.tar.gz test-suite

# Generating and uploading documentation

cd Doc

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

#
# Copy "Docs" to QE website

scp -R Doc <...>@<...>/wp-content/uploads/Doc-$version

# Connect to the website and create/update symbolic link to "Doc-$version"
