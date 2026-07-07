export QE='QE_Directory'

## polaron calculations
mpirun -np 32 ${QE}/pw.x < scf.in > scf.out
mpirun -np 32 ${QE}/pw.x < nscf.in > nscf.out
mpirun -np 32 ${QE}/epw.x -nk 32 < epw1.in > epw1.out
mpirun -np 1 ${QE}/epw.x -nk 1   < epw2.in > epw2.out
mpirun -np 32 ${QE}/epw.x -nk 32 < epw3.in > epw3.out

