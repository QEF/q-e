#!/bin/bash
home=`pwd`

cd example01/
./run_example
./compare.sh
echo "DONE"
rm -fr home/nicola/Scrivania/CODES/git/koopmans/quantum_espresso/qe_koopmans/tempdir
cd $home

echo "Running example02"
cd example02/
./run_example
./compare.sh
echo "DONE"
rm -fr home/nicola/Scrivania/CODES/git/koopmans/quantum_espresso/qe_koopmans/tempdir
cd $home 

echo "Running example03"
cd example03/
./run_example
./compare.sh
echo "DONE"
rm -fr home/nicola/Scrivania/CODES/git/koopmans/quantum_espresso/qe_koopmans/tempdir
cd $home 
