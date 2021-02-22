#!/bin/bash/

qepath="../../../../../../bin"
wannpath="../../../../../../wannier90-2.1.0/"

#sh clean.sh

mpirun -np 4 $qepath/pw.x < Si.scf.in > Si.scf.out 
mpirun -np 4 $qepath/pw.x < Si.nscf.in > Si.nscf.out 
echo " PWSCF DONE"


############# OCC STATE ###########
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


mpirun -np 4 $qepath/wann2kc.x < Si.kc_screen.in > Si.wann2kc.out
echo " Wann90 to KC DONE"

mpirun -np 4 $qepath/kc_screen.x < Si.kc_screen.in > Si.kc_screen.out
echo " SCREENING DONE"

printf "\n Screening parameters \n"
grep relaxed Si.kc_screen.out
echo " "

mpirun -np 4 $qepath/kc_ham.x < Si.kc_ham.in > Si.kc_ham.out
echo " HAMILTONIAN DONE"
