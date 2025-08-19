Author: Z. Dai

The example to calculate exciton polarons in LiF.

**For compiling QE and EPW: Please make sure that you turn on HDF5. To do so, usually you just need to do "./configure --with-hdf5=${YOUR_HDF5_DIRECTORY}", and then "make epw".
**For compiling BerkeleyGW: Also please make sure that you turn on HDF5. 

Below are the steps for running the calculations:

0. In setup.sh, specify the paths for QuantumEspresso and BerkeleyGW by setting $QEPATH and $BGWPATH.

1. Run the BSE calculations by running run.espresso.sh. This will calculate the BSE eigenvectors with zero exciton momentum. This step is to familiarize yourself with BSE calculations in BerkeleyGW.

2. In finiteQ_grid folder, run the BSE calculations for Q points in a uniform grid. 
   2.1 First run python3 generate_ex_band.py to prepare all the input files
   2.2 Then perform the BSE calculations by running run.bse.sh. 

3. In the epw folder, first run the scf and phonon calculations via run.ph.sh to get the variation of the scf potential.

4. Run pp.py to collect all the needed quantities.

5. IMPORTANT!!! Do NOT run another nscf calculations. Instead, just run run.epw.sh to perform the exciton polaron calculations.
   5.1 epw1.in is to generate the exciton-phonon matrix elements.
   5.2 epw2.in is to perform the exciton polaron calculations.
   5.3 epw3.in is to plot the electron and hole density of the exciton polaron.


Note: 
Sometimes in the finite Q calculations in BGW, the electronic wavefunctions normalize to a value very close to 1, but not 1. This is because BGW obtain the wave function at k+Q from the uniform k-grid as generated in the 05-Wfn_fi folder, but since what is imposed is the periodicity of \psi_nk in reciprocal space, the lattice periodic part u_nk is not actually periodic in reciprocal space. More specifically, u_nk+G_0 = e^i(G_0*r) u_nk. Therefore, if k+Q is beyond the first Brillouin zone, BGW will read the planewave coefficients of u_nk+Q-G_0 but map the G+G_0 component of it to the G component of u_nk+Q. Then, within a finite kinetic energy cutoff, some of the plane-wave components will be dropped at the boundary of the G-sphere because of the shift by G_0, and that’s way the wave functions do not normalize to 1 for some k+Q combinations.

The current example should run without a problem, but if you encounter such problems in other materials, try
(1) Increase the kinetic energy cut-off in the QE calculations and hope that one day when the cutoff is large enough, this problem will disappear 
(2) Hack into the BGW source code and change the threshold where they determine whether the wave functions normalize to 1. More specifically, in Common/misc.f90, search “if(abs(xnorm - 1.0d0) > TOL_Small)",  change TOL_Small to 1.0d-4, and recompile, and then the problems should disappear.
