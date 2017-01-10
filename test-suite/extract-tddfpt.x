#!/bin/sh
#
# Copyright (C) 2001-2016 Quantum ESPRESSO group
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# Maintainers: Filippo Spiga (filippo.spiga@quantum-espresso.org)
#              Samuel Ponce
#
# SP: This can lead to issue if you reach the OS pipe buffer 
#     You can increase the buffer in /proc/sys/fs/pipe-max-size  

  
fname=$1
args=$(echo $fname | awk -F= '{print $NF}')


### if [[ "$args" == "1" ]]
##  then
# SCF
e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`
### fi


# turbo_lanczos.x
alpha=`grep "alpha" $fname | awk '{print $2}'`
beta=`grep "beta " $fname | awk '{print $3}'`
gamma=`grep "gamma" $fname | awk '{print $2}'`


# turbo_spectrum.x
average=`grep "Average =" $fname | awk '{print $3}'`
averageosc=`grep "Average oscillation amplitude" $fname | awk '{print $5}'`
plotchi=`grep "chi_1_1=" *.plot_chi.dat | awk '{print $4}'`

if test "$plotchi" != ""; then
        echo plotchi
        for x in $plotchi; do echo $x; done
fi

if test "$average" != ""; then
        echo average
        for x in $average; do echo $x; done
fi

if test "$averageosc" != ""; then
        echo averageosc
        for x in $averageosc; do echo $x; done
fi

if test "$alpha" != ""; then
        echo alpha
        for x in $alpha; do echo $x; done
fi

if test "$beta" != ""; then
        echo beta
        for x in $beta; do echo $x; done
fi

if test "$gamma" != ""; then
        echo gamma
        for x in $gamma; do echo $x; done
fi

if test "$e1" != ""; then
	echo e1
	echo $e1
fi

if test "$n1" != ""; then
	echo n1
	echo $n1
fi

if test "$f1" != ""; then
	echo f1
	echo $f1
fi

if test "$p1" != ""; then
	echo p1
	echo $p1
