
# Example

Here we compute the energy current from one snapshot of SiO2 extracted from a CPMD simulation.
To calculate the energy current, for each snapshot at time t, we need to do the following:

- prepare the input, for example `input_energycurrent`:
```
 &energy_current
    delta_t=   0.500,
    init_linear= 'niente',
    file_output= 'current_hz',
    eta=   0.100,
    n_max=     5,
    status= 'compute',
 /
 &CONTROL
    calculation='md',
    restart_mode='from_scratch',
    pseudo_dir='./pseudo',
    outdir='./save',
    prefix= 'a_blocco',
    tprnfor=.true.,
/
 &SYSTEM
    ibrav=1,
    celldm(1)=   18.866000000000000,
    nat=    72,
    ntyp=     2,
    ecutwfc=  70.000,
/
 &ELECTRONS
    conv_thr = 1.D-10,
    mixing_beta = 0.7,
/
 & IONS
    ion_velocities = 'from_input',
/
 ATOMIC_SPECIES
   Si     28.08550000 Si_ONCV_PBE-1.1.upf
   O     15.99940000 O_ONCV_PBE-1.0.upf
 ATOMIC_POSITIONS {bohr}
   ...put here atomic positions, as in pw input file...

 ATOMIC_VELOCITIES
   ...put here atomic velocities...

 K_POINTS {Gamma}  
```
- run the code:
```
all_currents.x -in input_energycurrent
```

