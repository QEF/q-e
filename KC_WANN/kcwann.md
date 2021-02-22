project: KC_WANN
src_dir: ./src
output_dir: ./Doc_ford
summary: KC_WANN is an open-source F90/MPI code which calculates spectral properties using the koopmans-complinat framework and Maximally Localized Wannier Functions.
authors: Nicola Colonna
         Riccardo de Gennaro
         Nicola Marzari
predocmark: >
docmark_alt: #
predocmark_alt: <
display: public
         private
graph: true
source: false

KC_WANN stands for "Koopmans-compliant functionals in a Wannier representation". KC_WANN is an open-source F90/MPI 
code which calculates quasiparticle energies of finite and extended systems using
[Koopmans-compliant functionals](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.021051)
and [Maximally Localized Wannier Functions](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.56.12847).
The code consists of 3 programs:

1) an interface between Wannier90 and the KC_WANN code (wann_to_kc.x) 

2) a program that computes the screening coefficients as described [here](https://pubs.acs.org/doi/abs/10.1021/acs.jctc.7b01116) (kc_screen.x)

3) a program that computes the KC hamiltonian, eventually interpolate it and finally diagonalize it (kc_ham.x) 

KC_WANN is developed and maintained by [Nicola Colonna](https://people.epfl.ch/nicola.colonna) and [Riccardo de Gennaro](https://people.epfl.ch/riccardo.degennaro)

@Note
Still in development. This version of the code works with QE6.6

	TODO list: 
	1) Symmetry: at the moment no symmetry are used. 
	   This was primarely because to build up the periodic part of the Wannier
	   function [rho_of_q.f90] we need the KS wave-function on a regolar mesh
	   (pretty much like in W90). I see 2 possible strategy. The first one is to 
	   run PWscf with symmetry and then reconstruct all the KS wfcs on the regular
	   grid using the symmetry operations. The second is to find all the symmetry 
 	   and equivalent k/q points after the first part of the code. 
	2) Initialize the xc-kernel as the spin-polarized one also when nspin=1
	   This is needed to correctly define the perturbing potentials.
	3) A check on the WF spread to see if they are equivalent or not 
	   to speed up the calculation.

