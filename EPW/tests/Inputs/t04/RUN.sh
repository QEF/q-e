# Sequential run

../../../../bin/pw.x < scf.in > scf.out
../../../../bin/ph.x < ph.in > ph.out
python pp.py < pp.in
../../../../bin/pw.x < scf_epw.in > scf_epw.out
../../../../bin/pw.x < nscf_epw.in > nscf_epw.out
../../../src/epw.x < epw.in > epw.out

# Parallel run

mpirun -np 4 ../../../../bin/pw.x < scf.in > scf.out
mpirun -np 4 ../../../../bin/ph.x < ph.in > ph.out
python pp.py < pp.in
mpirun -np 4 ../../../../bin/pw.x < scf_epw.in > scf_epw.out
mpirun -np 4 ../../../../bin/pw.x -npool 4 < nscf_epw.in > nscf_epw.out
mpirun -np 4 ../../../src/epw.x -npool 4 < epw.in > epw.out
