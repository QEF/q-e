project: KCW
src_dir: ./src
output_dir: ./Doc_ford
summary: KCW is an open-source F90/MPI code which calculates spectral properties using the koopmans-compliant framework and Maximally Localized Wannier Functions.
authors: Nicola Colonna
         Riccardo de Gennaro
         Edward Linscott
         Nicola Marzari
predocmark: >
docmark_alt: #
predocmark_alt: <
display: public
         private
graph: true
source: false

#### NOTABENE To generate the automatic FORD documentation: "ford -I ../LAXlib kcw.md"

KCW stands for "Koopmans-compliant functionals in a Wannier representation". KCW is an open-source F90/MPI 
code which calculates quasiparticle energies of finite and extended systems using
[Koopmans-compliant functionals](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.021051)
and [Maximally Localized Wannier Functions](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.56.12847).
The details of this implementation are described [here](https://pubs.acs.org/doi/full/10.1021/acs.jctc.2c00161)
(see also [here](https://arxiv.org/abs/2202.08155) for the arXiv version). 
The code consists of 3 main modules specified by the "calculation" variable in the CONTROL namelist:

1) interface between PWscf, and Wannier90 and the KCW code (calculation="wann2kcw") 

2) calcuation of the screening coefficients (calculation="screen")

3) calculation, interpolation and diagonalization of the KC hamiltonian (calculation = "ham") 

Finally "calulation=cc" computes the (estimated) q+G=0 contribution to the bare and screened KC corrections.
A report on this quantities is printed on output and can be used to correct a posteriori a calculation performed 
without any corrective scheme (l_vcut=.false.)

KCW is developed and maintained by [Nicola Colonna](https://www.psi.ch/en/lns/people/nicola-colonna),  [Riccardo de Gennaro](https://people.epfl.ch/riccardo.degennaro), and [Edward Linscott](https://people.epfl.ch/edward.linscott)

	TODO list: 
	1) Symmetry: at the moment no symmetry are used. 
	This was primarely because W90 works on a regular mesh of k points and 
	does not account for symmetry. Prob need W90 to work in the IBZ.
	2) Initialize the xc-kernel and KC response as the spin-polarized one 
	also when nspin=1. This is needed to correctly define the perturbing
	potentials. At the moment we do a nspin=2 calculation from the beginning
	(not a big deal since the bottleneck is the LR calculation for which in 
        any case a spin-polarized response is needed).

