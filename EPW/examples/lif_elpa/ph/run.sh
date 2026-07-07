export QE='QE_Directory'

## phonon calculations
mpirun -np 32 ${QE}/pw.x < scf.in > scf.out
mpirun -np 32 ${QE}/ph.x < ph.in  > ph.out

