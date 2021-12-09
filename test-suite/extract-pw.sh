# Copyright (C) 2001 Quantum ESPRESSO
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

fname=$1

# maximum iterations extracted for relaxations
max_iter=3

# SCF
nks=`grep "number of Kohn-Sham states" $fname | awk '{print $5}'`
num_band=`awk "BEGIN{ print $nks * $max_iter }"`
e1=`grep ^! $fname | tail -1 | awk '{printf "%12.6f\n", $5}'`
n1=`grep 'convergence has' $fname | tail -1 | awk '{print $6}'`
f1=`grep "Total force" $fname | head -1 | awk '{printf "%8.4f\n", $4}'`
p1=`grep "P= " $fname | tail -1 | awk '{print $6}'`
band=`sed -n "/bands (ev)/{n;n;p}" $fname | awk '{print $1; print $2; print $3; print $4; print $5 }' | head -$num_band`

# NSCF
#ef1=`grep Fermi $fname | head -$max_iter | awk '{print $5}'`
ef1=$(awk 'BEGIN{ii=0} /^ *the Fermi energy is/{print $5; if(++ii>='$max_iter') exit;}' $fname)
eh1=`grep "highest occupied" $fname | tail -1 | awk '{print $5}'`
ehl1=`grep "highest occupied, lowest unoccupied" $fname | tail -1 | awk '{print $7; print $8}'`
tf1=`grep " P = " $fname | head -1 | awk '{printf "%7.5f", $3}'`

# extract geometry (volume excepted) after relaxation
geom=`sed -e '/new unit-cell/d' $fname |
  awk '/Begin final coordinates/,/End final coordinates/{
  # search all element in the line
  for ( i = 1; i<= NF; i++ ) {
    # print floating point numbers
    if(match($i, "[-+]?[0-9]+.?[0-9]+")!=0) {
      gsub(")", "", $i)
      print $i
    }
  }
}'`

vol=`grep 'new unit-cell' $fname | tail -1 | awk '{print $5}'`

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

if test "$ehl1" != ""; then
        echo ehl1
        for x in $ehl1; do echo $x; done
fi

if test "$band" != ""; then
        echo band
        for x in $band; do echo $x; done
fi

if test "$tf1" != ""; then
        echo tf1
        for x in $tf1; do echo $x; done
fi

if test "$vol" != ""; then
	echo vol
        echo $vol
fi

if test "$geom" != ""; then
	echo geom
        for x in $geom; do echo $x; done
fi
