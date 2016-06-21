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
  
fname=$1

# SCF
e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`

# NSCF
ef1=`grep Fermi $fname | awk '{print $5}'`
eh1=`grep "highest occupied" $fname | awk '{print $7}'`
el1=`grep "highest occupied" $fname | awk '{print $8}'`
tf1=`grep " P = " $fname | head -1 | awk '{printf "%7.5f", $3}'`

# EPW
q1=`grep "   q(" $fname | awk '{print $6 $7 $8}'`
dos1=`grep "DOS =" $fname | awk '{print $3}'`
e2=`grep " E(" $fname | awk '{print $4}'`
rsig=`grep "Re\[Sigma\]=" $fname | awk '{print $7}'` 
isig=`grep "Im\[Sigma\]=" $fname | awk '{print $10}'` 
z1=`grep " Z=" $fname | awk '{print $13}'`
lam=`grep "lam= " $fname | awk '{print $15}'`

if test "$q1" != ""; then
        echo q1
        for x in $q1; do echo $x; done
fi

if test "$dos1" != ""; then
        echo dos1
        echo $dos1
fi

if test "$e2" != ""; then
        echo e2
        for x in $e2; do echo $x; done
fi

if test "$rsig" != ""; then
        echo rsig
        for x in $rsig; do echo $x; done
fi

if test "$isig" != ""; then
        echo isig
        for x in $isig; do echo $x; done
fi

if test "$z1" != ""; then
        echo z1
        for x in $z1; do echo $x; done
fi

if test "$lam" != ""; then
        echo lam
        for x in $lam; do echo $x; done
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
fi


if test "$ef1" != ""; then
	echo ef1
	for x in $ef1; do echo $x; done
fi

if test "$eh1" != ""; then
        echo eh1
        for x in $eh1; do echo $x; done
fi

if test "$el1" != ""; then
        echo el1
        for x in $el1; do echo $x; done
fi

if test "$tf1" != ""; then
        echo tf1
        for x in $tf1; do echo $x; done
fi
