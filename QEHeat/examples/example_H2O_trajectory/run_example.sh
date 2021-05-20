ALL_CURRENTS=../../src/all_currents.x

mpirun -np 12 $ALL_CURRENTS -in input >& output
