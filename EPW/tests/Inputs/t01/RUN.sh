# Sequential run

../../../../bin/pw.x < scf.in > scf.out
../../../../bin/ph.x < ph.in > ph.out
python pp.py < pp.in
../../../../bin/pw.x < scf_epw.in > scf_epw.out
../../../../bin/pw.x < nscf_epw.in > nscf_epw.out
../../../src/epw.x < epw1.in > epw1.out
../../../src/epw.x < epw3.in > epw3.out
../../../src/epw.x < epw4.in > epw4.out
../../../src/epw.x < epw2.in > epw2.out

# Parallel run

mpirun -np 4 ../../../../bin/pw.x < scf.in > scf.out
mpirun -np 4 ../../../../bin/ph.x < ph.in > ph.out
python pp.py < pp.in
mpirun -np 4 ../../../../bin/pw.x < scf_epw.in > scf_epw.out
mpirun -np 4 ../../../../bin/pw.x -npool 4 < nscf_epw.in > nscf_epw.out
mpirun -np 4 ../../../src/epw.x -npool 4 < epw1.in > epw1.out
mpirun -np 4 ../../../src/epw.x -npool 4 < epw3.in > epw3.out
mpirun -np 4 ../../../src/epw.x -npool 4 < epw4.in > epw4.out
mpirun -np 4 ../../../src/epw.x -npool 4 < epw2.in > epw2.out
