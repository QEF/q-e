title: EPW
src_dir: ./src
output_dir: ./doc
project_website: http://epw.org.uk/
summary: EPW is the short name for "Electron-phonon Wannier". EPW is an open-source F90/MPI code which calculates properties related to the electron-phonon interaction using Density-Functional Perturbation Theory and Maximally Localized Wannier Functions.
authors: Samuel Poncé
         Roxana Margine
         Carla Verdi
         Feliciano Giustino
author_description: The EPW project is mainly developed at the university of Oxford.
github: https://github.com/sponce24
email: samuel.pon@gmail.com
project_sourceforge: http://qeforge.qe-forge.org/gf/project/q-e/
predocmark: >
media_dir: ./media
page_dir: ./Ford
docmark_alt: #
predocmark_alt: <
display: public
         private
source: false
graph: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: GNU
extra_filetypes: sh #

EPW is the short name for "Electron-phonon Wannier". EPW is an open-source F90/MPI
 code which calculates properties related to the electron-phonon interaction
 using [Density-Functional Perturbation Theory](http://journals.aps.org/rmp/abstract/10.1103/RevModPhys.73.515) 
and [Maximally Localized Wannier Functions](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.56.12847). 
EPW is developed and maintained by [Samuel Poncé](http://giustino.materials.ox.ac.uk/index.php/Site/SamuelPonc%e9), [Roxana Margine](http://www.binghamton.edu/physics/people/margine.html), [Carla Verdi](http://giustino.materials.ox.ac.uk/index.php/Site/CarlaVerdi), and [Feliciano Giustino](http://giustino.materials.ox.ac.uk/).

The reference technical manuscript for the latest pubic release is:
[EPW: Electron-phonon coupling, transport and superconducting properties using maximally localized Wannier functions](http://arxiv.org/abs/1604.03525)
by S. Poncé, E. R. Margine, C. Verdi, and F. Giustino. 


@Note
Since 26 April 2016 EPW is distributed as part of the [Quantum ESPRESSO](http://www.quantum-espresso.org/) suite. 

The code was written by Feliciano Giustino (EPW v1) and Jesse Noffsinger (EPW v2) while
 at the University of California, Berkeley. Brad Malone (Harvard) and Cheol-Hwan Park
 (Seoul National University) contributed with tests and benchmarks. 
Roxana Margine implemented the anisotropic Eliashberg theory while at the University of Oxford (EPW v3). 
Samuel Poncé (Oxford) made the code compatible with the latest version of Quantum Espresso v5 
in the latest release EPW v4. Carla Verdi (Oxford) developed the electron-phonon interpolation 
for polar materials including Froehlich correction (released within EPW v4).

EPW is based on the method introduced in F. Giustino et al, [Phys. Rev. B 76, 165108 (2007)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.76.165108). 
An extended description of the first public release has been published in 
J. Noffsinger et al, [Comput. Phys. Comm. 181, 2140 (2010)](http://www.sciencedirect.com/science/article/pii/S0010465510003218). The extension of EPW to include the
 anisotropic Midgal-Eliashberg theory is based on the method described in
 E. R. Margine et al, [Phys. Rev. B 87, 024505 (2013)](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.87.024505). 
The latest release of the code is described in S. Poncé et al, [arXiv:1604.03525](http://arxiv.org/abs/1604.03525). 
