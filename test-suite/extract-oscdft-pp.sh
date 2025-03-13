fname=$1
temp=`sed -n "/OSCDFT: eigenvalue/{n;p;}" $fname`
for e in $temp
do
    echo $e
done
