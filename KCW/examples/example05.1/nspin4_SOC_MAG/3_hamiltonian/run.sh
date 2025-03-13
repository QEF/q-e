#!/bin/bash
#para="mpirun -np 4"
printf "  KCW ham ... "
./link_wann.sh 
$para kcw.x -in kc.khi > kc.kho
echo "  done"

#printf "  Wann PostProc   ... "
#cd wannier_post 
#wannier90.x -pp wann
#pw2wannier90.x -in pw2wan.p2wi > pw2wan.p2wo
#wannier90.x wann
#echo "  done"
