export QE='QE_Directory'

## phonon calculations
#ibrun -np 32 ${QE}/pw.x < scf.in > scf.out
#ibrun -np 48 ${QE}/ph.x < ph.in  > ph.out

## dvscf postprocessing
## run pp.py


## epw calculations 
#ibrun -np 32 ${QE}/pw.x < nscf.in > nscf.out
#ibrun -np 32 ${QE}/epw.x -nk 32 < epw1.in > epw1.out

