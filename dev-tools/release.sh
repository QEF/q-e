#!/bin/sh -x

#tempdir=$HOME/Downloads
tempdir=/tmp
version=5.4.0

# make sure there is no locale setting creating unneeded differences.
LC_ALL=C
export LC_ALL

mkdir $tempdir
cd $tempdir
/bin/rm -rf espresso/ espresso-$version
# get the svn copy
svn checkout http://qeforge.qe-forge.org/svn/q-e/tags/QE-5.3.0/espresso
mv espresso/ espresso-$version/

cd espresso-$version

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
make doc

# generate PWGUI
make tar-gui PWGUI_VERSION=$version 
tar -xzvf PWgui-$version.tgz
/bin/rm PWgui-$version.tgz

cd ..

# core espresso

tar -czvf espresso-$version.tar.gz espresso-$version/archive \
                                   espresso-$version/clib \
                                   espresso-$version/configure \
                                   espresso-$version/COUPLE \
                                   espresso-$version/CPV \
                                   espresso-$version/dev-tools \
                                   espresso-$version/Doc \
                                   espresso-$version/environment_variables \
                                   espresso-$version/LAXlib \
                                   espresso-$version/FFTXlib \
                                   espresso-$version/Makefile \
                                   espresso-$version/include \
                                   espresso-$version/install \
                                   espresso-$version/License \
                                   espresso-$version/Modules \
                                   espresso-$version/PP \
                                   espresso-$version/pseudo \
                                   espresso-$version/PW \
                                   espresso-$version/README \
                                   espresso-$version/upftools


# Packages, ready for automatic unpacking

cd espresso-$version
tar -cvzf ../PWgui-$version.tar.gz    PWgui-$version
tar -czvf ../PHonon-$version.tar.gz   PHonon # PlotPhon QHA
tar -czvf ../EPW-$version.tar.gz      EPW
tar -czvf ../neb-$version.tar.gz      NEB
tar -czvf ../pwcond-$version.tar.gz   PWCOND
tar -czvf ../xspectra-$version.tar.gz XSpectra
tar -czvf ../GWW-$version.tar.gz      GWW
#tar -czvf ../GIPAW-$version.tar.gz    GIPAW
tar -czvf ../tddfpt-$version.tar.gz   TDDFPT
tar -czvf ../atomic-$version.tar.gz   atomic
tar -czvf ../test-suite-$version.tar.gz test-suite
tar -czvf ../EPW-$version.tar.gz EPW

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
