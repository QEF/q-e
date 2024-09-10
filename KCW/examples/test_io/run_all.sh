#!/bin/bash
home_local=`pwd`

cd example01/
./run_example
./compare.sh
echo "DONE"
cd $home_local

echo "Running example02"
cd example02/
./run_example
./compare.sh
echo "DONE"
cd $home_local 

echo "Running example03"
cd example03/
./run_example
./compare.sh
echo "DONE"
cd $home_local
