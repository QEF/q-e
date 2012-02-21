#!/bin/bash
###############################################
#This script moves cubes files without complaint so that the makefile is not broken
#OBM 2012
################################################


 if [ $# -ne 1 ] ; then
  echo "Usage:"
  echo "$0 <directory name>"
  exit
 fi
 echo "Moving current cube files to $1"
 if [ ! -d $1 ] ; then
  mkdir $1
 fi
 mv *.cube $1
 cd $1
 for file in *.cube
 do
  gzip $file
 done
 cd ..

 
