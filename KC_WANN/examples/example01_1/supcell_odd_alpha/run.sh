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

nbnd=`grep num_bands Si.win | awk '{print $3}'`
nbnd_sqr=`echo $nbnd | awk '{ print $1*$1}'`
echo $nbnd > nbnd
tail -$nbnd_sqr Si_u.mat > mat
cat nbnd mat > unimatrx_occ_up.txt
rm  mat


############# EMPT STATE ###########
$wannpath/wannier90.x -pp Si_emp.win
mpirun -np 4 $qepath/pw2wannier90.x < Si_emp.pw2wann.in > Si_emp.pw2wann.out
echo " PW2WANN EMP DONE"

$wannpath/wannier90.x  Si_emp.win
echo " WANN EMP DONE"

nbnd_wann=`grep num_wann Si_emp.win | awk '{print $3}'`
nbnd_ks=`grep num_bands Si_emp.win | awk '{print $3}'`
nbnd_sqr=`echo $nbnd_wann $nbnd_ks | awk '{ print $1*$2}'`
echo $nbnd_wann $nbnd_ks > nbnd
tail -$nbnd_sqr Si_emp_u.mat > mat
cat nbnd mat > unimatrx_empt_up.txt
rm  mat

