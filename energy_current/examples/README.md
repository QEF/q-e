###
############ INSTALLATION OVERVIEW #############
###
### Before running the examples, you need to be sure that QEHeat is installed.
### To do this, first perform a standard Quantum-Espresso installation. For a basic one, in the main directory execute the commands:
### ./configure
### make pw 
### This should produce the standard pw.x executable. To compile Qeheat just you need to enter the correct folder :
### cd /energy_current/src
### make
### This should produce the executable all_currents.x in the src and bin folders.
###
############ EXAMPLES #########################
### 
### See also ../../Doc/INPUT_ALL_CURRENTS.html for a description of the inputs.
### To run an example, enter in the respective trajectory and execute the script run_example.sh. Modify it if necessary.
### In general, running a Qeheat calculation just needs the execution of the command: all_currents.x -in input_energycurrent , 
### after the input file has been prepared.
###
### Each example comes with a reference folder where the output files can be compared with the ones produced by a new installation/run.
###
### Example 1 and 2 need a parallel installation to finish in a reasonable time. Example 1 was run in the reference calculation on 4 cores and Example 2 on 96.
### Example 3 can be easily run on a single core implementation.



# 
# 1. example_SiO2_single  
# 



Here we compute the energy current and its indivial components for a single snapshot of Silica. 
For this purpose, one needs the additional namelist in the input file :

&energy_current
    delta_t=   0.500,
    file_output= 'current_hz',
    eta=   0.100,
    n_max=     5,
 /

& IONS
    ion_velocities = 'from_input',

and the CARD ATOMIC_VELOCITIES in the input file must be filled as well with the istantenous atomic velocities.
These are the ingredients needed for a basic single snapshot calculation.
 
The files produced, apart from the standard output, are :

- file_output.dat : this reports the total energy current, the eletronic current and the center of mass velocity for each species.
Only file_output.dat needs to be used to evaluate the thermal conductivity coefficient.

- file_output : this reports the total energy current divided in individual components, if a more specific analysis is needed.

Each file has a header specifying the output units. See also the description ../../Doc/INPUT_ALL_CURRENTS.html 



# 2. example_H2O_trajectory



Here we compute the energy current from a previously computed Car-Parrinello trajectory which is provided together with the input file . The trajectory provided
here comes from a 125 water molecule simulation.

To calculate the energy current, for every timestep of the trajectory located in  ${trajdir}.pos` and ${trajdir.vel}, (velocities are in cp units in this example) 
we need to insert some additional keywords in the energy_current namelists :

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

note that the only differences are `trajdir` and `first_step` in the `energy_current` namelist. Stil the IONS namelist needs to be set 
ion_velocities='from_input' and the ATOMIC_VELOCITIES card must be filled. Nevertheless, note that if first_step=1 the snapshot of the input file will be skipped.
Check also the last_step keyword in the documentation.

The output with the istantenous energy currents is written in the file specified by `file_output` and `file_output.dat`, as in example 1. 
The same format of the single snapshot calculation is kept, with data from all the snapshots of the trajectory appended sequantially.



# 3. example_small_H2O_trajectory



This example is very similar to the previous one, but the trajectory is computed on-the-fly via the cp.x program of the installed distribution. It produces the  
trajectory of a single water molecule and therefore the calculation is suited for a serial environment.




