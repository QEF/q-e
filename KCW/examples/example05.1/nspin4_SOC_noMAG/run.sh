#!/bin/bash

echo "DFT calculation"
cd 0_dft
 ./clean.sh
 ./run.sh
cd ../

echo "KI calculation"
echo " 1-Wannier"
cd 1_wannier
 ./clean.sh
 ./run.sh
cd ../

echo " 2-Screening"
cd 2_screening
./clean.sh
./run.sh
cd ..

echo " 3-Hamiltonian"
cd 3_hamiltonian
./clean.sh
./run.sh
cd ..

