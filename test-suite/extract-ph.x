# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
 
fname=$1

# SCF
e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`

# PH
diel=`grep -A 4 '  Dielectric constant in cartesian' $fname | grep -v '  Dielectric constant' | awk '{print $2; print $3; print $4 }'`
born=`grep "     E[x-z]  ( " $fname | awk '{print $3; print $4; print $5}'`
phfreq=`grep "     freq (.*THz" $fname | awk '{print $5; print $8}'`



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

if test "$diel" != ""; then
        echo diel
        for x in $diel; do echo $x; done
fi
 
if test "$born" != ""; then
        echo born
        for x in $born; do echo $x; done
fi

if test "$phfreq" != ""; then
        echo phfreq
        for x in $phfreq; do echo $x; done
fi


