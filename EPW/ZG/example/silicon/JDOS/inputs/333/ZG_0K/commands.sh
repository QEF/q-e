awk '/eigenvalues size/{x=NR+44}(NR<=x){print}' silicon.xml > oo
sed -i '/eigenvalues size/d' oo # remove lines containin a string
sed -i '/^$/d' oo # remve empty lines
tr '\n' ' ' < oo > kk
tr -s ' '  '\n' < kk > data_one_column.dat
sed -i '/^$/d' data_one_column.dat # remve empty lines
rm oo kk
