#!/bin/bash/

qepath="../../../../bin"
wannpath="../../../../wannier90-2.1.0/"

#sh clean.sh

mpirun -np 4 $qepath/pw.x < Si.scf.in > Si.scf.out 
echo " PWSCF DONE"


$wannpath/wannier90.x -pp Si.win 
mpirun -np 4  $qepath/pw2wannier90.x < Si.pw2wann.in > Si.pw2wann.out 
echo " PW2WANN OCC DONE"

$wannpath/wannier90.x Si.win 
echo " WANN OCC DONE"

############# EMPT STATE ###########
$wannpath/wannier90.x -pp Si_emp.win
mpirun -np 4 $qepath/pw2wannier90.x < Si_emp.pw2wann.in > Si_emp.pw2wann.out
echo " PW2WANN EMP DONE"

$wannpath/wannier90.x  Si_emp.win
echo " WANN EMP DONE"

mpirun -np 4 $qepath/wann_to_kc.x < Si.kc_screen_occ.in > Si.wann_to_kc.out
echo " Wann90 to KC OCC DONE"

#mpirun -np 4 $qepath/kc_screen.x < Si.kc_screen_occ.in > Si.kc_screen_occ.out
#echo " ALPHA OCC DONE"

#printf "\n Relevant info \n"
#grep relaxed Si.kc_screen_occ.out
#echo " "

#mpirun -np 4 $qepath/kc_screen.x < Si.kc_screen_emp.in > Si.kc_screen_emp.out
#echo " ALPHA EMP DONE"
#
#printf "\n Relevant info \n"
#grep relaxed Si.kc_screen_emp.out
#echo " "



