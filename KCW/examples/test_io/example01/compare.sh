#!/bin/bash
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

tol_ene=0.000001
tol_eig=0.001
tol_r2=0.00001
tol_a=0.0001
check=1
echo 
echo " Testing the results ..." 
echo "  tolerances: "
printf "     energy %10.7f Ry  \n" $tol_ene
printf "     eigval %10.7f eV  \n" $tol_eig
printf "     R^2    %10.7f A^2 \n" $tol_r2
printf "     alpha  %10.7f     \n" $tol_a

for iors in ".true." ".false."; do
for iosp in  ".true." ".false."; do

echo "  Testing the different IO schema against ../../example01/reference"
echo "  Testing io_real_space = $iors io_single_precision = $iosp" 

## CHECK alpha parameters
check=1
for i in `seq 1 8`; do
    a=`grep relaxed results/iors_${iors}_iosp_$iosp/Si.kcw-screen.out | head -$i | tail -1|  awk '{print $12}'`
a_ref=`grep relaxed ../../example01/reference/Si.kcw-screen.out | head -$i | tail -1 | awk '{print $12}'`
err=`echo $i $a $a_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
if (( $(echo "$err > $tol_a" |bc -l) )); then
 echo -e "  ${RED}KC_SCREEN WARNING: ${NC}alpha does not match for orbital $i: $a $a_ref $err"
 check=0
fi
done
if (( $check )); then echo -e "  ${GREEN}KC_SCREEN   OK!${NC}"; fi

check=1

for i in `seq 1 8`; do 
 grep "KI  " results/iors_${iors}_iosp_$iosp/Si.kcw-ham.out | head -$i | tail -1 > pp
 grep "KI  " ../../example01/reference/Si.kcw-ham.out | head -$i | tail -1 >> pp
 for j in `seq 1 8`; do 
  col=`echo $j+1 | bc`
  eig=`head -1 pp | awk -v col=$col '{print $col}'`
  eig_ref=`tail -1 pp | awk -v col=$col '{print $col}'`
  err=`echo $i $eig $eig_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
  if (( $(echo "$err > $tol_eig" |bc -l) )); then
   echo -e "  ${RED}KC_HAM WARNING: ${NC}eig does not match ik=$i,v=$j: $eig $eig_ref $err"
   check=0
  fi
 done
done
if (( $check )); then echo -e "  ${GREEN}KC_HAM      OK!${NC}"; fi
rm pp
echo 

done
done
