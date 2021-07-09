#!/bin/bash
declare -i i j
file1="si.333_ZG_bands" #without .in
file2="kpoints.dat"
#
i=1
wc -l "$file2" > aa
a=$(awk '{print $1}' aa)
rm aa
echo "Kpts are" $a
#
while [ $i -le $a ];do
  $(awk 'NR == '$i' ' $file2  > testmz) # print correct k-point 
  cat $file1 testmz > "$file1"_"$i".in
  i=$((i+1))
done
rm testmz
