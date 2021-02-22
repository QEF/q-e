#!/bin/bash/
rm -fr *out

for charge in 0.0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1; do 
occ=`echo 1.0-$charge | bc -l`
echo $charge $occ

cat > h2o.scf_charge$charge.in << EOF
&CONTROL
  calculation='scf',
  restart_mode='from_scratch',
  pseudo_dir = '../../pseudo'
  prefix = 'h2o'
  outdir ='./out'
  wf_collect=.true.
/
&SYSTEM
  ecutwfc =   45.0
  ibrav = 0
  input_dft = 'PBE'
  nat = 3
  nspin = 2
  ntyp = 2
  nbnd = 4
  assume_isolated='mt'
  occupations='from_input'
  tot_charge = $charge
/
&ELECTRONS
    diagonalization='david'
    mixing_mode = 'plain'
    mixing_beta = 0.7
    conv_thr =  0.5d-12
    startingpot = 'file'
/
ATOMIC_SPECIES
H 1 H_ONCV_PBE-1.0.upf 
 O 1 O_ONCV_PBE-1.0.upf 

ATOMIC_POSITIONS angstrom
 O 6.7571 6.0000 5.9023166667 
 H 7.5142 6.0000 6.4884166667 
 H 6.0000 6.0000 6.4884166667

CELL_PARAMETERS angstrom
9.5142 0.0 0.0 
0.0 8.0 0.0 
0.0 0.0 8.5861 
K_POINTS automatic
1 1 1 0 0 0
OCCUPATIONS
1.0 1.0 1.0 $occ
1.0 1.0 1.0 1.0
EOF
mpirun -np 4 ../../../../bin/pw.x < h2o.scf_charge$charge.in > h2o.scf_charge$charge.out 
echo " charge $charge DONE"

done


