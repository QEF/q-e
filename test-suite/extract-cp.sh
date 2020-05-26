# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
 
fname=$1
e1=`grep "total energy =" $fname | tail -1 | awk '{printf "%18.6f\n", $4}'`
s1=`grep -A 3 "Total stress" $fname | tail -3 | tr '\n' ' ' | awk '{ printf "%-18.8f", $1+$2+$3+$4+$5+$6+$7+$8+$9 }'`
v1u=`grep -A 2 "Eigenvalues (eV).*spin.*1" $fname | tail -1 | awk '{ for(i=1;i<=NF;i++) { v=v+$i; } print v }'` 
v1d=`grep -A 2 "Eigenvalues (eV).*spin.*2" $fname | tail -1 | awk '{ for(i=1;i<=NF;i++) { v=v+$i; } print v }'` 
t1=`grep -A 6 "Averaged Physical Quantities"  $fname | tail -1 | awk '{ print $4 }'`
 
if test "$e1" != ""; then
	echo e1
	echo $e1
fi

if test "$s1" != ""; then
	echo s1
	echo $s1
fi

if test "$v1u" != ""; then
	echo v1u
	echo $v1u
fi

if test "$v1d" != ""; then
	echo v1d
	echo $v1d
fi

if test "$t1" != ""; then
	echo t1
	echo $t1
fi
