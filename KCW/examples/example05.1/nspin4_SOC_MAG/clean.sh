#!/bin/bash

rm -fr out 
cd 0_dft
 ./clean.sh
cd ../

cd 1_wannier
 ./clean.sh
cd ../

cd 2_screening
./clean.sh
cd ..

cd 3_hamiltonian
./clean.sh
cd ..

