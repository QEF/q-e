#!/bin/bash

python3 rotate.py inpos rpos invel rvel
python3 rotate.py rmat 

cat > all_currents.in << EOF
&energy_current
    delta_t=   0.500,
    file_output= 'current_hz',
    eta=   0.100,
    n_max=     5,
 /
 &CONTROL
    calculation='md',
    restart_mode='from_scratch',
    pseudo_dir='../pseudo',
    outdir='./save',
    disk_io = 'none',
/
 &SYSTEM
    ibrav=0,
    celldm(1) = 10.,
    nat=    3,
    ntyp=     2,
    ecutwfc=  80.0,
/
 &ELECTRONS
    conv_thr = 1.D-11,
    mixing_beta = 0.7,
/
 &IONS
    ion_velocities = 'from_input',
/
 ATOMIC_SPECIES
   H      1.00000000 H_HSCV_PBE-1.0.upf
   O     16.00000000 O_HSCV_PBE-1.0.upf
ATOMIC_POSITIONS {bohr}
$(paste -d '' types rpos)

ATOMIC_VELOCITIES
$(paste -d '' types rvel)

CELL_PARAMETERS {alat}
$(cat rmat)

K_POINTS {Gamma}
EOF
rm rpos rvel rmat


cat > all_currents_not_rotated.in << EOF
&energy_current
    delta_t=   0.500,
    file_output= 'current_hz_not_rotated',
    eta=   0.100,
    n_max=     5,
 /
 &CONTROL
    calculation='md',
    restart_mode='from_scratch',
    pseudo_dir='../pseudo',
    outdir='./save',
    disk_io = 'none',
/
 &SYSTEM
    ibrav=1,
    celldm(1) = 10.,
    nat=    3,
    ntyp=     2,
    ecutwfc=  80.0,
/
 &ELECTRONS
    conv_thr = 1.D-11,
    mixing_beta = 0.7,
/
 &IONS
    ion_velocities = 'from_input',
/
 ATOMIC_SPECIES
   H      1.00000000 H_HSCV_PBE-1.0.upf
   O     16.00000000 O_HSCV_PBE-1.0.upf
ATOMIC_POSITIONS {bohr}
$(paste -d '' types inpos)

ATOMIC_VELOCITIES
$(paste -d '' types invel)

K_POINTS {Gamma}
EOF

