# Sequential run

../../../../bin/pw.x < scf.in > scf.out
../../../../bin/ph.x < ph.in > ph.out
python pp.py < pp.in
../../../../bin/pw.x < scf_epw.in > scf_epw.out
../../../../bin/pw.x < nscf_epw.in > nscf_epw.out
../../../src/epw.x < epw_iso.in > epw_iso.out
../../../src/epw.x < epw_iso_real.in > epw_iso_real.out
../../../src/epw.x < epw_aniso.in > epw_aniso.out

# Parallel run

mpirun -np 4 ../../../../bin/pw.x < scf.in > scf.out
mpirun -np 4 ../../../../bin/ph.x < ph.in > ph.out
python pp.py < pp.in
mpirun -np 4 ../../../../bin/pw.x < scf_epw.in > scf_epw.out
mpirun -np 4 ../../../../bin/pw.x -npool 4 < nscf_epw.in > nscf_epw.out
mpirun -np 4 ../../../src/epw.x -npool 4 < epw_iso.in > epw_iso.out
mpirun -np 4 ../../../src/epw.x -npool 4 < epw_iso_real.in > epw_iso_real.out
mpirun -np 4 ../../../src/epw.x -npool 4 < epw_aniso.in > epw_aniso.out
