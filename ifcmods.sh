#!/bin/sh

echo "creating catalog file for Intel compiler v.6 and earlier..."

topdir=`pwd`

pcl_ph="D3"
pcl_pw="PH PP Gamma PWCOND pwtools $pcl_ph"
pcl_modules="PW CPV flib upftools atomic $pcl_pw"
pcl_dot="Modules $pcl_modules"

for dir in $pcl_dot
do
  echo work.pc > $topdir/$dir/intel.pcl
done

for dir in $pcl_modules
do
   echo ../Modules/work.pc >> $topdir/$dir/intel.pcl
done

for dir in $pcl_pw
do
   echo ../PW/work.pc >> $topdir/$dir/intel.pcl
done

for dir in $pcl_ph
do
   echo ../PH/work.pc >> $topdir/$dir/intel.pcl
done
echo "done"

