#export QE='QE_Directory'

## phonon calculations
ibrun -np 32 ${QE}/pw.x < scf.in > scf.out
ibrun -np 32 ${QE}/ph.x < ph.in  > ph.out

