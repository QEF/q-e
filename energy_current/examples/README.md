
# Example single

Here we compute the energy current from one snapshot.
To calculate the energy current, for each snapshot at time t, we need to do the following:

- prepare the input, for example `input_energycurrent`:
```
 &energy_current
    delta_t=   0.500,
    file_output= 'current_hz',
    eta=   0.100,
    n_max=     5,
 /
 &CONTROL
    calculation='md',
    restart_mode='from_scratch',
    pseudo_dir='./pseudo',
    outdir='./save',
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
   ...put here atomic velocities... (usual PW units)

 K_POINTS {Gamma}  
```
- run the code:
```
all_currents.x -in input_energycurrent
```

# Example trajectory

Here we compute the energy current from a cp trajectory
To calculate the energy current, for every timestep of the trajectory `traj.pos` and `traj.vel` (velocities are in cp units in this example) located in the folder 'trajectory', we need to do the following:

- prepare the input, for example `input_energycurrent`:
```
 &energy_current
    delta_t=   0.500,
    file_output= 'current_hz',
    eta=   0.100,
    n_max=     5,
    trajdir='trajectory/traj'
    first_step=1,
    vel_input_units='CP'
 /
 &CONTROL
    calculation='md',
    restart_mode='from_scratch',
    pseudo_dir='./pseudo',
    outdir='./save',
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
   ...put here atomic velocities... (CP units)

 K_POINTS {Gamma}  
```
- run the code:
```
all_currents.x -in input_energycurrent
```

note that the only differences are `trajdir` and `first_step` in the `energy_current` namelist. The snapshot of the input file will be skipped.

# Output

The output is written in the file specified by `file_output`. Then a column formatted file is written in `file_output`.dat, with the following content:
```
STEP t_ps Jx Jy Jz
.    .    .  .  .
.    .    .  .  .
.    .    .  .  .
```
STEP and t_ps are the same of the input trajectory.
