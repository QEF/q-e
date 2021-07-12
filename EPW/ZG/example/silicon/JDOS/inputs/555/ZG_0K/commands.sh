array_mz2=( 1 2 3 4 ) #100 200 300 400 500 600 750 1000 30 50 10 75 ) #6 7 9 14
for j in ${array_mz2[@]}
   do
   awk '/eigenvalues size/{x=NR+150}(NR<=x){print}' silicon_"$j".xml > oo
   sed -i '/eigenvalues size/d' oo # remove lines containin a string
   sed -i '/^$/d' oo # remve empty lines
   tr '\n' ' ' < oo > kk
   tr -s ' '  '\n' < kk > data_one_column_"$j".dat
   sed -i '/^$/d' data_one_column_"$j".dat # remve empty lines
   rm oo kk
done
