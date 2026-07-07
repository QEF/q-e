Author: K. Luo
Date: Jan 26, 2026

> A minimal test case to compare electron-phonon coupling (g-matrix) elements 
  computed by PH and EPW for LiF.

### Below are the steps for running the calculations:
1. In ph_g/, run the SCF and PH calculations by executing run.sh. 
2. Run SCF and PH calculations in epw/, this generates the dvscf files for the electron-phonon matrix elements.
3. Run pp.py in EPW/bin/ to collect all required quantities.
4. In epw/, run SCF, NSCF, and EPW1 calculations by exexcuting run.sh.
5. Compare the g-matrix elements of selected band indices (m,n) and phonon mode (nu). 

### Note:
1. In both ph.out and epw1.out, one may see multiple entries for a g-matrix element:
   "|g_sym|[meV]   |g|[meV]   Re(g)[meV]   Im(g)[meV]". i
   Here, the first entry "|g_sym|" is the symmetrized magitude obtained by algebraically 
   averaging raw magnitude |g| over degenerate electronic states (m,n) and phonon mode (nu).
   The real and imaginary parts reported correspond to the unsymmetrized complex g 
   and are included for completeness.
2. In EPW output, the g-matrix elements include an acoustic sum rule (ASR) correction 
   applied to the phonon eigenmodes. As a result, values in ph.out and epw1.out can differ 
   slightly, especially for acoustic modes near Gamma.
   By switching 'asr_typ' to 'no' in epw/epw1.in, these discrepancies can be closed.
3. The complex phase of g depends on both electronic wavefunction gauge and phonon eigenmodes. 
   Since these can differ between PH and EPW, the complex g values are generally not directly 
   comparable; comparisons should focus on gauge-invariant quantities |g_sym|.


