#!/bin/bash
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

tol_ene=0.000001
tol_eig=0.001
tol_r2=0.00001
tol_sh=0.00001
check=1
echo 
echo " Testing the results ..." 
echo "  tolerances: "
printf "     energy       %10.7f Ry  \n" $tol_ene
printf "     eigval       %10.7f eV  \n" $tol_eig
printf "     R^2          %10.7f A^2 \n" $tol_r2
printf "     Self-Hartree %10.7f Ry \n" $tol_sh


for iors in ".true." ".false."; do
for iosp in  ".true." ".false."; do

echo "  Testing the different IO schema against ../../example03/reference"
echo "  Testing io_real_space = $iors io_single_precision = $iosp"

check=1

for i in `seq 1 8`; do 
 grep "KI  " results/iors_${iors}_iosp_${iosp}/Si.kcw-ham.out | head -$i | tail -1 > pp
 grep "KI  " ../../example03/reference/Si.kcw-ham.out | head -$i | tail -1 >> pp
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
if (( $check )); then echo -e "  ${GREEN}KC_HAM        OK!${NC}"; fi
rm pp

check=1

for i in `seq 1 8`; do 
 grep "SH  " results/iors_${iors}_iosp_${iosp}/Si.kcwpp_sh.out | head -$i | tail -1 > pp
 grep "SH  " ../../example03/reference/Si.kcwpp_sh.out | head -$i | tail -1 >> pp
 for j in `seq 1 8`; do 
  col=`echo $j+1 | bc`
  sh=`head -1 pp | awk -v col=$col '{print $col}'`
  sh_ref=`tail -1 pp | awk -v col=$col '{print $col}'`
  err=`echo $i $sh $sh_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
  if (( $(echo "$err > $tol_sh" |bc -l) )); then
   echo -e "  ${RED}KC_HAM WARNING: ${NC}eig does not match ik=$i,v=$j: $sh $sh_ref $err"
   check=0
  fi
 done
done
if (( $check )); then echo -e "  ${GREEN}KCWPP_SH      OK!${NC}"; fi
rm pp

check=1

for i in `seq 1 62`; do 
 grep -A 2 "KC interpolate" results/iors_${iors}_iosp_${iosp}/Si.kcwpp_interp.out | grep -v "KC inter" | grep -v "\-\-" | sed -r '/^\s*$/d' | head -$1 | tail -1 > pp
 grep -A 2 "KC interpolate" ../../example03/reference/Si.kcwpp_interp.out | grep -v "KC inter" | grep -v "\-\-" | sed -r '/^\s*$/d' | head -$1 | tail -1 > pp
 for j in `seq 1 8`; do 
  col=`echo $j+1 | bc`
  eig=`head -1 pp | awk -v col=$col '{print $col}'`
  eig_ref=`tail -1 pp | awk -v col=$col '{print $col}'`
  err=`echo $i $eig $eig_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
  if (( $(echo "$err > $tol_eig" |bc -l) )); then
   echo -e "  ${RED}KCWPP_INTERP WARNING: ${NC}eig does not match ik=$i,v=$j: $eig $eig_ref $err"
   check=0
  fi
 done
done
if (( $check )); then echo -e "  ${GREEN}KCWPP_INTERP  OK!${NC}"; fi
rm pp
echo

done
done
