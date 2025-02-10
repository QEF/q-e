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
# SCF'
e1=`grep ^! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`
### fi

## PP WANNIER
nkp1=`grep "     1     8   " $fname | tail -1 | awk '{print $3}'`
nkp2=`grep "     1     8   " $fname | tail -1 | awk '{print $4}'`
nkp3=`grep "     1     8   " $fname | tail -1 | awk '{print $5}'`

## PW2WANN 
eig1=`head -1  $fname | awk '{print $3}'`

## WANNIER 
omega=`grep "Omega Total" $fname | awk '{print $7}'`

## WANN2KCW
sh=`grep "orb     1" $fname | awk '{print $4}' | tail -1`

##KCW screen
rpi=`grep "iwann  =     1" $fname | awk '{print $6}'`
upi=`grep "iwann  =     1" $fname | awk '{print $9}'`

##KCW hamilt
homo_ks=`grep "KS       highest occupied level (ev):" $fname | tail -1 | awk '{print $6}'`
homo_kc=`grep "KI\[2nd\]  highest occupied level (ev):" $fname | tail -1 | awk '{print $6}'`
gap_ks=`grep "KS       highest occupied, lowest unoccupied level (ev):" $fname | tail -1 | awk '{print $9-$8}'`
gap_kc=`grep "KI\[2nd\]  highest occupied, lowest unoccupied level (ev):" $fname | tail -1 | awk '{print $9-$8}'`

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

if test "$nkp1" != ""; then
	echo nkp1
	echo $nkp1
fi

if test "$nkp2" != ""; then
	echo nkp2
	echo $nkp2
fi

if test "$nkp3" != ""; then
	echo nkp3
	echo $nkp3
fi

if test "$eig1" != ""; then
	echo eig1
	echo $eig1
fi
if test "$omega" != ""; then
	echo omega
	echo $omega
fi
if test "$sh" != ""; then
	echo sh
	echo $sh
fi
if test "$rpi" != ""; then
	echo rpi
	echo $rpi
fi
if test "$upi" != ""; then
	echo upi
	echo $upi
fi
if test "$homo_ks" != ""; then
	echo homo_ks
	echo $homo_ks
fi
if test "$homo_kc" != ""; then
	echo homo_kc
	echo $homo_kc
fi
if test "$gap_ks" != ""; then
	echo gap_ks
	echo $gap_ks
fi
if test "$gap_kc" != ""; then
	echo gap_kc
	echo $gap_kc
fi
