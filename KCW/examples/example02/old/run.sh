#!/bin/bash/
rm -fr *out
mpirun -np 4 ../../../bin/pw.x < h2o.scf.in > h2o.scf.out 
echo " PWSCF DONE"
mpirun -np 4 ../../../bin/wann2kc.x < h2o.kc_screen.in > h2o.wann2kc.out 
echo " Wann90 to KC DONE"
mpirun -np 4 ../../../bin/kc_screen.x < h2o.kc_screen.in > h2o.kc_screen.out 
echo " KC screen DONE"

printf "\n Relevant info \n"
grep relaxed h2o.kc_screen.out
echo " "

mpirun -np 4 ../../../bin/kc_ham.x < h2o.kc_ham.in > h2o.kc_ham.out
echo " KC ham DONE"

