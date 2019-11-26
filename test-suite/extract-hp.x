# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.
#
# SP: This can lead to issue if you reach the OS pipe buffer
#     You can increase the buffer in /proc/sys/fs/pipe-max-size
  
fname=$1
args=$(echo $fname | awk -F= '{print $NF}')

##echo $fname > /home/sponce/program/espresso/test-suite/tmp.txt 
###    echo $args >> /home/sponce/program/espresso/test-suite/tmp.txt

### if [[ "$args" == "1" ]]
##  then
# SCF
e1=`grep ! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`
### fi


# HP
u=`grep "Hubbard U (eV)" -A 1 $fname | tail -n 1 | awk '{print $7}'`


if test "$u" != ""; then
        echo u
        echo $u
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

