#!/bin/bash/

qepath="../../../../../../bin"
wannpath="../../../../../../bin/"

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


mpirun -np 4 $qepath/kcw.x < Si.kcw-wann2kcw.in > Si.kcw-wann2kcw.out
echo " Wann90 to KC DONE"

mpirun -np 4 $qepath/kcw.x < Si.kcw-screen.in > Si.kcw-screen.out
echo " SCREENING DONE"

printf "\n Screening parameters \n"
grep relaxed Si.kcw-screen.out
echo " "

mpirun -np 4 $qepath/kcw.x < Si.kcw-ham.in > Si.kcw-ham.out
echo " HAMILTONIAN DONE"
