
This example illustrates how to use pw.x and ph.x to calculate phonon
frequencies at Gamma and X for Si and C in the diamond structure and 
for fcc-Ni.

The calculation proceeds as follows (for the meaning of the cited input
variables see the appropriate INPUT_* file)

1) make a self-consistent calculation (input=si.scf.in,
   output=si.scf.out)

2) make a phonon calculation for gamma point (input=si.phG.in,
   output=si.phG.out). Note that outdir and prefix have the same
   values as in the input for scf calculation.
   Note that the q-point (the Gamma point in this case) is read
   after the namelist inputph. 
   Dynamical matrix will be saved for further analysis in fildyn (whose
   default name is 'matdyn', here overwritten to si.dynG)

   The Gamma point is the ONLY one for which the phonon code can be run
   directly after the scf run (provided save file was defined).
   For all other points in the BZ a non-scf preparatory run is needed,
   as it is shown here for the X point.

3) make a phonon calculation for X point (input=si.phX.in,
   output=si.phX.out). Note that the q-point (the X point 
   in this case) is read after the namelist &inputph.
   Dynamical matrices will be saved for further analysis in fildyn
   (whose default name is 'matdyn', here overwritten to si.dynX). 
   The fildyn file contains the dynamical matrices for q-points that
   are equivalent by symmetry to the q-point given in input (in this
   case the three X points), therefore when calculating phonons in the
   whole BZ the code needs to be run only on inequivalent points.

   The phonon code can perform also a 'single-mode' calculation where
   only a single atom is moved in a single direction and the corresponding
   column of the dynamical matrix is calculated. This may be useful when
   studying vibrational properties of defects and/or surfaces.
   This feature is illustrated in the following part of the test.

4) Perform again a scf calculation for Si. Why? because the phonon
   calculation performed in the previous step has overwritten the
   original wavefunctions. This step can be avoided by saving the
   content of outdir and restoring it after step 3)

5) make a phonon calculation at X point for a single mode
   (input=si.phXsingle.in, output=si.phXsingle.out). 
   Apart from the specification of the mode to be calculated
   (modenum=3), the input is identical to case 4). The only difference
   is that fildyn is not given a name, so the default one (matdyn) is
   used. This is done in order not to overwrite the previous files. 
   The total dynamical matrix is not computed and matdyn is actually 
   left empty.

6) make a scf calculation of C (described with US-PP) in the diamond
   structure. 
   (input=c.scf.in, output=c.scf.out).

7) make a phonon calculation at Gamma point for C in the diamond structure. 
   Note that epsil=.true. and the response to an electric field and 
   the dielectric constant are also calculated.
   (input=c.phG.in, output=c.phG.out).

8) make a scf calculation at of magnetic fcc-Ni.
   (input=ni.scf.in, output=ni.scf.out).

9) make a phonon calculations at X point for magnetic fcc-Ni.
   (input=ni.phX.in, output=ni.phX.out). 
