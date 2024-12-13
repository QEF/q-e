#!/bin/bash
home=`pwd`

cd example01/
./run_example
./compare.sh
echo "DONE"
cd $home

echo "Running example01_sym"
cd example01_sym/
./run_example
./compare.sh
echo "DONE"
cd $home

echo "Running example02"
cd example02/
./run_example
./compare.sh
echo "DONE"
cd $home 

echo "Running example03"
cd example03/
./run_example
./compare.sh
echo "DONE"
cd $home 

echo "Running example04"
cd example04/
./run_example
./compare.sh
echo "DONE"
cd $home 

echo "Running example05"
cd example05/
./run_example
./compare.sh
echo "DONE"
cd $home

echo "Running example06"
cd example06/
./run_example
./compare.sh
echo "DONE"
cd $home


