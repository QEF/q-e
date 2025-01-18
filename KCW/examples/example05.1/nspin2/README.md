# Description

This example show how to use pw.x and kcw.x to compute the 
electronic structure of FCC Silicon. As a reference the band
structure at the LDA level is first computed in 0_dft.
The entire KI calculation is split down into three steps:

* 1_init: Contains the inputs for the initial DFT calculation and for the interface with the KCW code
* 2_screening: Contains the input for the calculation of the KI screening coefficients (on a coarse k-mesh)
* 3_hamiltonian: Contains the inputs for the final KI calculation
