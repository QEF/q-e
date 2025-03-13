#!/bin/bash

if [ -e "all.f90" ] 
then
    echo "removing all.f90"
    rm all.f90
fi

cat *f90 > all.f90
make  # need mod files
cp Makefile Makefile.orig
cp Makefile.all Makefile
make
mv Makefile.orig Makefile

# pick out external subs
nm all.o | grep -v " T " | grep -v " B " | grep -v "_gfortran" | grep -v "MOD" | grep "_" | grep -v " b " > ext_deps
# don't need to highlight lapack/blas calls
cat ext_deps | grep -vi "zgemm"  | grep -vi "zhpevx"  | grep -vi "zgemv"  | grep -vi "zdotc" | grep -vi "zaxpy" | grep -vi "ddot" > ext_deps1
mv ext_deps1 ext_deps
cat ext_deps