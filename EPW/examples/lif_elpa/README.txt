Author: K. Luo
Date: Jan 06, 2026

> The example to calculate the visualize excited-state polarons in LiF.

### For compiling ELPA:
Please make sure ELPA is properly compiled before compiling QE/EPW if you wish to accelerate large-scale calculations. Flags and directories may be machine-dependent. For example, one can configure ELPA in Permutter as:
"./configure --enable-gpu -prefix=/pscratch/sd/s/sabyadk/QE_with_ELPA/elpa/elpa/ 'FCFLAGS=-O2' 'CFLAGS=-O2' 'FFLAGS= -allow' FC=ftn CC=cc CXX=CC -disable-avx512 --disable-openmp --enable-static --with-GPU-compute-capability=sm_60 --enable-c-tests=no "

### For compiling QE and EPW:
To integrate ELPA solvers, configure QE/EPW as:
"./configure --enable-openmp --enable-parallel --with-scalapack=intel --with-elpa-lib=$ELPA_LIBS/libelpa_openmp.a --with-elpa-include=$ELPA_INCLUDE "

Then build EPW: "make epw "

Make sure the ELPA directories are set correctly, e.g.:
"export ELPA_LIBS='.../elpa/lib' "
"export ELPA_INCLUDE='.../elpa/include/elpa_openmp-2021.11.002/modules' "

### Below are the steps for running the calculations:
1. In the ph/ folder, run SCF and phonon calculations by executing run.sh to obtain the variation of the SCF potential.
2. Run pp.py in EPW/bin/ to collect all required quantities.
3. In the epw/ folder, run SCF and NSCF calculations:
   - 3.1 epw1.in: generate the electron-phonon matrix elements.
   - 3.2 epw2.in: perform ground-state polaron calculations. A number of excited states are also calculated at each iteration; the number of polaron eigenstates is controlled by 'nstate_plrn'.
   - 3.3 epw3.in: plot the (excited-state) polaron charge density. The state index is controlled by 'plot_psir_plrn'. The sign is retained if 'lsign_psir_plrn' is turned on.

### Note:
1. If the ELPA solver is not integrated, 'eigen_solver_plrn = "elpa"' cannot be used.
2. The relaxation of an excited-state polaron in epw2.in can be tested by changing 'istate_relax_plrn' (use at your own risk).
3. Setting 'prtvkk = .true.' in epw2.in stores velocity matrix elements on the fine grid, which is necessary if oscillator strengths of polarons are needed.
4. In the epw/ folder, run clean.sh to remove output files and restart the EPW workflow.

