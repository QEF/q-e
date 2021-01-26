![q-e-logo](logo.jpg)

This is the distribution of the Quantum ESPRESSO suite of codes (ESPRESSO:
opEn-Source Package for Research in Electronic Structure, Simulation, and
Optimization)

[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

## USAGE

Quick installation instructions for the impatient, using "make":
(`[]` means "optional"):
```
./configure [options]
make all
```
"make" alone prints a list of acceptable targets. Optionally,
`make -jN` runs parallel compilation on `N` processors.
Link to binaries are found in bin/.

Using "CMake":

```
mkdir ./build
cd ./build
cmake [-DCMAKE_INSTALL_PREFIX=/path/to/install] ..
make [-jN]
[make install]
```
"make" builds all targets. Link to binaries are found in build/bin.
If `make install` is invoked, directory `CMAKE_INSTALL_PREFIX`
is prepended onto all install directories.

For more information, see the general documentation in directory Doc/, 
package-specific documentation in \*/Doc/, and the web site 
http://www.quantum-espresso.org/. Documentation for developers 
can be found on [Wiki page on gitlab](https://gitlab.com/QEF/q-e/-/wikis/home).

## PACKAGES

- PWscf: structural optimisation and molecular dynamics on the electronic ground state, with self-consistent solution of DFT equations;
- CP: Car-Parrinello molecular dynamics;
- PHonon: vibrational and dielectric properties from DFPT (Density-Functional Perturbation Theory);
- TD-DFPT: spectra from Time-dependent DFPT;
- HP: calculation of Hubbard parameters from DFPT;
- EPW: calculation of electron-phonon coefficients in metals;
- PWCOND: ballistic transport;
- XSpectra: calculation of X-ray absorption spectra;
- PWneb: reaction pathways and transition states with the Nudged Elastic Band method;
- GWL: many-body perturbation theory in the GW approach using ultra-localised Wannier functions and Lanczos chains.

## Modular libraries
The following libraries have been isolated and partially encapsulated in view of their release for usage in other codes as well:

- UtilXlib: performing basic MPI handling, error handling, timing handling.
- FFTXlib: parallel (MPI and OpenMP) distributed three-dimensional FFTs, performing also load-balanced distribution of data (plane waves, G-vectors and real-space grids) across processors.
- LAXlib: parallel distributed dense-matrix diagonalization, using ELPA, SCALapack, or a custom algorithm.
- KS Solver: parallel iterative diagonalization for the Kohn-Sham Hamiltonian (represented as an operator),using block Davidson and band-by-band or block Conjugate-Gradient algorithms.
- LRlib: performs a variety of tasks connected with (time-dependent) DFPT, to be used also in connection with Many-Body Perturbation Theory.
- upflib: pseudopotential-related code

## GPU-enabled version
Since Feb.2021 this repository also works for GPU's (currently only NVidia). See file [README_GPU.md](README_GPU.md).

## Contributing
Quantum ESPRESSO is an open project: contributions are welcome.
Read the [Contribution Guidelines](CONTRIBUTING.md) to see how you
can contribute.

## LICENSE

All the material included in this distribution is free software;
you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation;
either version 2 of the License, or (at your option) any later version.

These programs are distributed in the hope that they will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
675 Mass Ave, Cambridge, MA 02139, USA.
