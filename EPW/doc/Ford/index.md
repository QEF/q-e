title: EPW overview
author: Samuel Poncé
date: 01-06-2016


<div style="text-align:center"><img src ="http://epw.org.uk/figures/logo_v7.png" width="600"></div>


EPW is an open-source F90/MPI code which calculates properties related to the electron-phonon interaction
 using [Density-Functional Perturbation Theory](http://journals.aps.org/rmp/abstract/10.1103/RevModPhys.73.515)
and [Maximally Localized Wannier Functions](http://journals.aps.org/prb/abstract/10.1103/PhysRevB.56.12847).

EPW is licensed under a [[GNU General Public License]](http://www.gnu.org/licenses/gpl-3.0.en.html)

EPW is part of the [[Quantum Espresso]](http://www.quantum-espresso.org/) software package.

## Installation 

The EPW software is only tested and intened to run on Linux (Mac OS might work but not tested).

* Download the latest version of [[Quantum-ESPRESSO]](http://www.qe-forge.org/gf/project/q-e/frs/?action=FrsReleaseBrowse&frs_package_id=18).

* Unpack and configure Quantum-ESPRESSO
```bash
tar -xvf espresso-5.4.0.tar.gz && cd espresso-5.4.0 && ./configure
``` 

* Compile EPW (this will also compile pwscf, phonon, and wannier90)
```bash
make -j 4 pwall
make -j 4 ph
make -j 4 epw
```

* The executable will be available in espresso-5.4.0/bin/epw.x or espresso-5.4.0/EPW/bin/epw.x
