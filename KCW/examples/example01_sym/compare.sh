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

## CHECK PWSCF
a=`grep ! results/Si.scf.out | awk '{print $5}'`
a_ref=`grep ! reference/Si.scf.out | awk '{print $5}'`
err=`echo $a $a_ref | awk '{printf "%20.15f \n", sqrt(($1-$2)*($1-$2))}'`
if (( $(echo "$err > $tol_ene" |bc -l) )); then
 echo -e "${RED}  PWSCF  WARNING: ${NC} Total energy does not match: $a $a_ref $err"
 check=0
fi
if (( $check )); then echo -e "  ${GREEN}PWSCF       OK!${NC}"; fi

check=1
a=`grep high results/Si.scf.out | awk '{print $7}'`
a_ref=`grep high reference/Si.scf.out | awk '{print $7}'`
err=`echo $a $a_ref | awk '{printf "%20.15f \n", sqrt(($1-$2)*($1-$2))}'`
if (( $(echo "$err > $tol_eig" |bc -l) )); then
 echo -e "${RED}  PWNSCF WARNING: ${NC} VBM does not match: $a $a_ref $err"
 check=0
fi
a=`grep high results/Si.scf.out | awk '{print $8}'`
a_ref=`grep high reference/Si.scf.out | awk '{print $8}'`
err=`echo $a $a_ref | awk '{printf "%20.15f \n", sqrt(($1-$2)*($1-$2))}'`
if (( $(echo "$err > $tol_eig" |bc -l) )); then
 echo -e "${RED}  PWNSCF WARNING: ${NC} CBM does not match: $a $a_ref $err"
 check=0
fi
if (( $check )); then echo -e "  ${GREEN}PWNSCF      OK!${NC}"; fi



check=1

## CHECK W90 R^2
##   OCC
for i in `seq 2 5`; do
    a=`grep -A 4 "Final State" results/wann_occ/Si.wout | head -$i | tail -1|  awk '{print $11}'`
a_ref=`grep -A 4 "Final State" reference/wann_occ/Si.wout | head -$i | tail -1 | awk '{print $11}'`
#echo $i $a $a_ref | awk '{print $1-1, $2, $3, $2-$3}'
iorb=`echo $i | awk '{print $1-1}'`
err=`echo $i $a $a_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
if (( $(echo "$err > $tol_r2" |bc -l) )); then
# echo "WARNING: R^2 does not match for orbital $iorb: $a $a_ref $err"
 echo -e "${RED}  W90 WARNING: ${NC} R^2 does not match for orbital $iorb: $a $a_ref $err"
 check=0
fi
done
##   EMP
for i in `seq 2 5`; do
    a=`grep -A 4 "Final State" results/wann_emp/Si_emp.wout | head -$i | tail -1|  awk '{print $11}'`
a_ref=`grep -A 4 "Final State" reference/wann_emp/Si_emp.wout | head -$i | tail -1 | awk '{print $11}'`
#echo $i $a $a_ref | awk '{print $1-1+4, $2, $3, $2-$3}'
iorb=`echo $i | awk '{print $1-1+4}'`
err=`echo $i $a $a_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
if (( $(echo "$err > $tol_r2" |bc -l) )); then
 echo -e "${RED}  W90 WARNING: ${NC} R^2 does not match for orbital $iorb: $a $a_ref $err"
 check=0
fi
done
if (( $check )); then echo -e "  ${GREEN}WANNIER     OK!${NC}"; fi

## CHECK alpha parameters
check=1
for i in `seq 1 8`; do
    a=`grep relaxed results/Si.kcw-screen.out | head -$i | tail -1|  awk '{print $12}'`
a_ref=`grep relaxed reference/Si.kcw-screen.out | head -$i | tail -1 | awk '{print $12}'`
err=`echo $i $a $a_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
if (( $(echo "$err > $tol_a" |bc -l) )); then
 echo -e "  ${RED}KC_SCREEN WARNING: ${NC}alpha does not match for orbital $i: $a $a_ref $err"
 check=0
fi
done
if (( $check )); then echo -e "  ${GREEN}KC_SCREEN   OK!${NC}"; fi

check=1

for i in `seq 1 8`; do 
 grep "KI  " results/Si.kcw-ham.out | head -$i | tail -1 > pp
 grep "KI  " reference/Si.kcw-ham.out | head -$i | tail -1 >> pp
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
