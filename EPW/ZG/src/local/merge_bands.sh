#!/bin/bash
declare -i i j
file1="bands01.dat" #without .in
file2="kpoints.dat"
file3="spectral_weights01.dat"
i=1
ls kpt_*/bands01.dat > oo
wc -l oo > aa
a=$(awk '{print $1}' aa)
rm aa
mkdir all_bands
cd all_bands
touch aa1
touch aa2
#
echo "Kpts are" $a
echo "&plot nbnd= 135, nks=     "$a" /" > aa1
echo "&plot nbnd= 135, nks=     "$a" /" > aa2
#
    while [ $i -le $a  ];do
        cp ../kpt_$i/$file1 "$file1"_"$i"
        sed -i '1,1d' "$file1"_"$i"
        cat aa1 "$file1"_"$i"  > oo
        mv oo aa1
        rm "$file1"_"$i"
# for spectral weights
        cp ../kpt_$i/$file3 "$file3"_"$i"
        sed -i '1,1d' "$file3"_"$i"
        cat aa2 "$file3"_"$i"  > oo
        mv oo aa2
        rm "$file3"_"$i"
        i=$((i+1))
    done
mv aa1 $file1
mv aa2 $file3
