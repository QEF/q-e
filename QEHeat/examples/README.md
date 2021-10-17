# INSTALLATION OVERVIEW

 Before running the examples, you need to be sure that QEHeat is installed.
 To compile, we suggest the following procedure. As usual in a QuantumEspresso installation, enter the distribution folder and run autoconf :

```
 > ./configure
```

 Than from the distribution folder:

```
 >  make all_currents
```

 This will produce the executable `all_currents.x`, the executable for QEHeat, in the respective src and bin folders.


# EXAMPLES

 See also `../Doc/INPUT_ALL_CURRENTS.html` for a description of the inputs.
 To run an example, enter in the respective trajectory and execute the script `run_example.sh`. Modify it if necessary.
 In general, running a Qeheat calculation just needs the execution of the command: 
  `all_currents.x -in input_energycurrent` , 
 after the input file has been prepared.

 Each example comes with a reference folder where the output files can be compared with the ones produced by a new installation/run.
 Pseudopotentials can be downloaded from [http://www.quantum-simulation.org/potentials/sg15_oncv/]()
 For the examples, the following pseudos were used:
 H_HSCV_PBE-1.0.upf,  O_HSCV_PBE-1.0.upf,  O_ONCV_PBE-1.0.upf,  Si_ONCV_PBE-1.1.upf 
 , which should be present in the folder `pseudo`

 Example 1 and 2 need a parallel installation to finish in a reasonable time. Example 1 was run in the reference calculation on 4 cores and Example 2 on 12. 

 Example 3 can be easily run on a single core (serial) installation. Note that Example 3 requires the program `cp.x` to be installed, for example with `make cp`
 from the main distribution folder.



# 1. example_SiO2_single  



Here we compute the energy current and its indivial components for a single snapshot of Silica. 
For this purpose, one needs the additional namelist in the input file :

```
&energy_current
    delta_t=   0.500,
    file_output= 'current_hz',
    eta=   0.100,
    n_max=     5,
 /

& IONS
    ion_velocities = 'from_input',
```

and the CARD ATOMIC_VELOCITIES in the input file must be filled as well with the istantenous atomic velocities.
These are the ingredients needed for a basic single snapshot calculation.
 
The files produced, apart from the standard output, are :

- `file_output.dat` : this reports the total energy current, the eletronic current and the center of mass velocity for each species.
Only file_output.dat needs to be used to evaluate the thermal conductivity coefficient.

- `file_output` : this reports the total energy current divided in individual components, if a more specific analysis is needed. This file is mainly thought for development purposes.

`file_output.dat` comes with a header specifying the output units. The same units are used in the more detailed current decomposion given in `file_output`.
 See also the description ../Doc/INPUT_ALL_CURRENTS.html 



# 2. example_H2O_trajectory



Here we evaluate the energy current from a previously computed Car-Parrinello (CP) trajectory, which is provided together with the input file. The trajectory provided 
comes from a 125 water molecule simulation.

We calculate the energy current for every timestep of the trajectory located in  `${trajdir}.pos` and `${trajdir}.vel` (velocities are in CP units in this example). 
For this purpose we need to insert some additional keywords in the energy_current namelists :

```
 &energy_current
    delta_t=   0.500,
    file_output= 'current_hz',
    eta=   0.100,
    n_max=     5,
    trajdir='traj'
    first_step=1,
    vel_input_units='CP'
 /
```

note that the only different keywords with respect to a single snapshot calculations are `trajdir` and `first_step` in the `energy_current` namelist. Still, in the IONS namelist 
the keyword ion_velocities='from_input' must be set and the ATOMIC_VELOCITIES card must be filled.

The output with the istantenous energy currents is written in the files with names `${file_output}` and `${file_output}.dat`, as in example 1. 
The same format of the single snapshot calculation is kept, data from all the snapshots of the trajectory being appended sequentially.

If computational time is an issue, for this example we suggest to insert the flag last_step=xxx in the energy_current namelist,
replacing xxx with the desired index step, and than execute the run script.
Note that the keywords  `last_step` and  `first_step` refer to the indexes reported in the files  `${trajdir}.pos` and `${trajdir}.vel` and are not sequential indexes. The snapshot
in the input file is assigned an index 0. As a concrete example, the combination of `first_step=1`  and `last_step=953008` will skip the snapshot of the input file and
evaluate only the first snapshot of the trajectory because 953008 is the first index that appears in the trajectory file.



# 3. example_small_H2O_trajectory



This example is very similar to the previous one, but a Car-Parrinello trajectory is computed on-the-fly via the cp.x program of the just installed QE distribution. It produces the  
trajectory of a single water molecule and therefore the calculation is suited for a serial environment. 

This example requires the program cp.x to be installed. If this is not the case, you can enter the distribution folder and run "make cp".

Note that the trajectory produced by cp.x will be probably different due to the stochasticity inherent in the Car-Parrinello molecular dynamics simulation. For exact comparison 
with a novel installation one can substitute `trajdir='reference/traj/cp'` and comment in the run_example_water script the call to cp.x. This way the files produced by `file_output`
should be comparable with the reference, up to numerical noise that is always present in the finite difference derivative with non perfectly converged wavefunctions.




