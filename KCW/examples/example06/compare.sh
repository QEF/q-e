#!/bin/bash
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

tol_ene=0.000001
tol_eig=0.001
tol_r2=0.00001
tol_a=0.0001

neigtotest=10
nstart=`echo 42-$neigtotest | bc`

check=1
echo 
echo " Testing the results ..." 
echo " Testing the $neigtotest highest KCW eigenvalues ..." 
echo "  tolerances: "
printf "     energy %10.7f Ry  \n" $tol_ene
printf "     eigval %10.7f eV  \n" $tol_eig
printf "     R^2    %10.7f A^2 \n" $tol_r2
printf "     alpha  %10.7f     \n" $tol_a

## CHECK PWSCF
a=`grep ! results/cri3.scf.out | awk '{print $5}'`
a_ref=`grep ! reference/cri3.scf.out | awk '{print $5}'`
err=`echo $a $a_ref | awk '{printf "%20.15f \n", sqrt(($1-$2)*($1-$2))}'`
if (( $(echo "$err > $tol_ene" |bc -l) )); then
 echo -e "${RED}  PWSCF  WARNING: ${NC} Total energy does not match: $a $a_ref $err"
 check=0
fi
if (( $check )); then echo -e "  ${GREEN}PWSCF (Etot)     OK!${NC}"; fi

check=1
a=`grep high results/cri3.scf.out | awk '{print $7}'`
a_ref=`grep high reference/cri3.scf.out | awk '{print $7}'`
err=`echo $a $a_ref | awk '{printf "%20.15f \n", sqrt(($1-$2)*($1-$2))}'`
if (( $(echo "$err > $tol_eig" |bc -l) )); then
 echo -e "${RED}  PWNSCF WARNING: ${NC} VBM does not match: $a $a_ref $err"
 check=0
fi
a=`grep high results/cri3.scf.out | awk '{print $8}'`
a_ref=`grep high reference/cri3.scf.out | awk '{print $8}'`
err=`echo $a $a_ref | awk '{printf "%20.15f \n", sqrt(($1-$2)*($1-$2))}'`
if (( $(echo "$err > $tol_eig" |bc -l) )); then
 echo -e "${RED}  PWNSCF WARNING: ${NC} CBM does not match: $a $a_ref $err"
 check=0
fi
if (( $check )); then echo -e "  ${GREEN}PWSCF (eig)      OK!${NC}"; fi


echo "  SPIN-UP ..." 

## CHECK W90 R^2
for ib in 1 2 3 4 5; do
check=1
if [ $ib = 1 ]; then nw=2;
elif [ $ib = 2 ]; then nw=6;
elif [ $ib = 3 ]; then nw=6
elif [ $ib = 4 ]; then nw=24
elif [ $ib = 5 ]; then nw=4
fi
grep -A $nw "Final State" results/SPIN_UP/wann_block$ib/wann.wout    | tail -n $nw | awk '{print $11}' | sort -n > spread_res
grep -A $nw "Final State" reference/SPIN_UP/wann_block$ib/wann.wout    | tail -n $nw | awk '{print $11}' | sort -n > spread_ref
for i in `seq 1 $nw`; do
    a=`cat spread_res | head -$i | tail -1|  awk '{print $11}'`
a_ref=`cat spread_ref | head -$i | tail -1 | awk '{print $11}'`
#echo $i $a $a_ref | awk '{print $1-1, $2, $3, $2-$3}'
iorb=`echo $i | awk '{print $1-1}'`
err=`echo $i $a $a_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
if (( $(echo "$err > $tol_r2" |bc -l) )); then
# echo "WARNING: R^2 does not match for orbital $iorb: $a $a_ref $err"
 echo -e "  ${RED}WANNIER   BLOCK $ib  WARNING: ${NC} R^2 does not match for orbital $iorb: $a $a_ref $err"
 check=0
fi
done
if (( $check )); then echo -e "    ${GREEN}WANNIER   BLOCK $ib  OK!${NC}"; fi
done
rm spread_res spread_ref


## CHECK SH kcw interface

check=1 
grep "    SH " results/SPIN_UP/cri3.kcw-wann2kcw.out | awk '{print $4, $2}' | sort -n > sh_res
grep "    SH " reference/SPIN_UP/cri3.kcw-wann2kcw.out | awk '{print $4, $2}' | sort -n > sh_ref
for i in `seq 1 42`; do
     a=`cat sh_res | head -$i | tail -1 | awk '{print $4*13.6057}'`
 a_ref=`cat sh_ref | head -$i | tail -1 | awk '{print $4*13.6057}'`
iorb=$i 
err=`echo $i $a $a_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
if (( $(echo "$err > $tol_eig" |bc -l) )); then
 echo -e "  ${RED}KCW-w2kcw  WARNING: ${NC} Self-Hartree does not match for orbital $iorb: $a $a_ref $err"
 check=0
fi
done 
if (( $check )); then echo -e "    ${GREEN}KCW-w2kcw          OK!${NC}"; fi
rm sh_res sh_ref




check=1

#printf "  ik "
for i in `seq 1 8`; do 
# printf "$i "
 ii=`echo $i*6 | bc`
 grep "KI  " results/SPIN_UP/cri3.kcw-ham.out | head -$ii | tail -6 | sed -e "s/KI//" | xargs > pp
 grep "KI  " reference/SPIN_UP/cri3.kcw-ham.out | head -$ii | tail -6 | sed -e "s/KI//" | xargs >> pp
 for j in `seq $nstart 42`; do 
  col=$j
  eig=`head -1 pp | awk -v col=$col '{print $col}'`
  eig_ref=`tail -1 pp | awk -v col=$col '{print $col}'`
  err=`echo $i $eig $eig_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
  if (( $(echo "$err > $tol_eig" |bc -l) )); then
   echo -e "  ${RED}KCW-ham WARNING: ${NC}eig does not match ik=$i,v=$j: $eig $eig_ref $err"
   check=0
  fi
 done
done
if (( $check )); then echo -e "    ${GREEN}KCW-ham            OK!${NC} (tested $neigtotest higest eigenvalues)"; fi
rm pp

check=1



echo "  SPIN-DOWN ..." 

## CHECK W90 R^2
for ib in 1 2 3 4 5 6; do
check=1
if [ $ib = 1 ]; then nw=2;
elif [ $ib = 2 ]; then nw=6;
elif [ $ib = 3 ]; then nw=6
elif [ $ib = 4 ]; then nw=18
elif [ $ib = 5 ]; then nw=6
elif [ $ib = 6 ]; then nw=4
fi
grep -A $nw "Final State" results/SPIN_DOWN/wann_block$ib/wann.wout    | tail -n $nw | awk '{print $11}' | sort -n > spread_res
grep -A $nw "Final State" reference/SPIN_DOWN/wann_block$ib/wann.wout    | tail -n $nw | awk '{print $11}' | sort -n > spread_ref
for i in `seq 1 $nw`; do
    a=`cat spread_res | head -$i | tail -1|  awk '{print $11}'`
a_ref=`cat spread_ref | head -$i | tail -1 | awk '{print $11}'`
#echo $i $a $a_ref | awk '{print $1-1, $2, $3, $2-$3}'
iorb=`echo $i | awk '{print $1-1}'`
err=`echo $i $a $a_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
if (( $(echo "$err > $tol_r2" |bc -l) )); then
# echo "WARNING: R^2 does not match for orbital $iorb: $a $a_ref $err"
 echo -e "  ${RED}WANNIER   BLOCK $ib  WARNING: ${NC} R^2 does not match for orbital $iorb: $a $a_ref $err"
 check=0
fi
done
if (( $check )); then echo -e "    ${GREEN}WANNIER   BLOCK $ib  OK!${NC}"; fi
done
rm spread_res spread_ref

## CHECK SH kcw interface

check=1
grep "    SH " results/SPIN_DOWN/cri3.kcw-wann2kcw.out | awk '{print $4, $2}' | sort -n > sh_res
grep "    SH " reference/SPIN_DOWN/cri3.kcw-wann2kcw.out | awk '{print $4, $2}' | sort -n > sh_ref
for i in `seq 1 42`; do
     a=`cat sh_res | head -$i | tail -1 | awk '{print $4*13.6057}'`
 a_ref=`cat sh_ref | head -$i | tail -1 | awk '{print $4*13.6057}'`
iorb=$i
err=`echo $i $a $a_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
if (( $(echo "$err > $tol_eig" |bc -l) )); then
 echo -e "  ${RED}KCW-w2kcw  WARNING: ${NC} Self-Hartree does not match for orbital $iorb: $a $a_ref $err"
 check=0
fi
done
if (( $check )); then echo -e "    ${GREEN}KCW-w2kcw          OK!${NC}"; fi
rm sh_res sh_ref


check=1 
#printf "  ik "
for i in `seq 1 8`; do
# printf "$i "
 ii=`echo $i*6 | bc`
 grep "KI  " results/SPIN_DOWN/cri3.kcw-ham.out | head -$ii | tail -6 | sed -e "s/KI//" | xargs > pp
 grep "KI  " reference/SPIN_DOWN/cri3.kcw-ham.out | head -$ii | tail -6 | sed -e "s/KI//" | xargs >> pp
 for j in `seq $nstart 42`; do
  col=$j
  eig=`head -1 pp | awk -v col=$col '{print $col}'`
  eig_ref=`tail -1 pp | awk -v col=$col '{print $col}'`
  err=`echo $i $eig $eig_ref | awk '{printf "%20.15f \n", sqrt(($2-$3)*($2-$3))}'`
  if (( $(echo "$err > $tol_eig" |bc -l) )); then
   echo -e "  ${RED}KC-ham WARNING: ${NC}eig does not match ik=$i,v=$j: $eig $eig_ref $err"
   check=0
  fi
 done
done
if (( $check )); then echo -e "    ${GREEN}KCW-ham            OK!${NC} (tested $neigtotest higest eigenvalues)"; fi
rm pp
echo
