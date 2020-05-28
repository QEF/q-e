ALL_CURRENTS=../../src/all_currents.x

mpirun -np 4 $ALL_CURRENTS -in input_energycurrent >& output_energycurrent
