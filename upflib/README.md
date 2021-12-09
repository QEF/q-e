
# Library of pseudopotential code

This directory contains a library of pseudopotential-related code,
extracted from the Quantum ESPRESSO distribution. This library depends only
upon module mp.f90 of UtilXlib and upon a few modules and routines of devXlib;
upon a few LAPACK routines; requires a suitable `../make.inc` file in Makefile. 
Other than this, it can be independently compiled.

Currently, it includes
- basic definitions of the UPF (Unified Pseudopotential File) format
- basic I/O operations on UPF files
- setup of the interpolation tables and of other basic variables 
- interpolation of pseudopotentials
- generation of various pseudopotentials matrix elements
- utilities: spherical harmonics and Bessel functions, integration routines
Old UPF specifications can be found here:
http://www.quantum-espresso.org/pseudopotentials/unified-pseudopotential-format
The xml schema for the newer UPF definition can be found here:
http://www.quantum-espresso.org/ns/qes/qe_pp-1.0.xsd

In addition to the `libupf.a` library, executable utilities are produced:

- `upfconv.x`, converting pseudopotentials in other formats into UPF:
   see `upfconv.x -h` for more

- `virtual_v2.x`, courtesy Jingyang Wang (jw598@cornell.edu), generates
   an averaged pseudopotential suitable for Virtual Crystal Approximation

- `casino2upf.x`, courtesy Mike Towler (see below)

A python script `fixfile.py` is also present, to remove undesired `&`
characters from UPF files that hinder their parsing by xml tools.

## CASINO and QE pseudopotentials

The following notes are kept for reference (they might be obsolete).
Code `upfconv.x -c` should replace code `upf2casino2.x` mentioned below.
Code `casino2upf.x` was moved to upflib/ and works (?) again, at least
for the example provided by Jake Muff, since v.6.8. Old notes start here:

Two utilities are provided with the Quantum Espresso distribution to 
enable the PWscf code to be used in conjunction with the CASINO quantum 
Monte Carlo code.

Of course all pseudopotentials generated via these automatic tools should 
be tested before being used for production runs.

It should be noted that ultrasoft and PAW pseudopotentials cannot be used
with the CASINO code. Currently only UPF files containing norm-conserving 
pseudopotentials can be converted using these utilities.

### casino2upf.x

The first of these is casino2upf.x . This utility takes a given CASINO 
tabulated pseudopotential file and one or more awfn.data files specifying 
the pseudoatomic wavefunctions to be used in creating the 
Kleinman-Bylander projectors. A UPF file containing the projectors and the 
local potential is then written to the file name specified in inputpp. Any
errors are communicated to the user via stderr.

Usage:
	
        ./casino2upf.x < inputpp

A sample inputpp file for converting a Trail and Needs pseudopotential 
would be:

```
inputpp:
	&inputpp
		pp_data='pp.data'
		upf_file='my_pseudo_potential.UPF'
	/
	3
	awfn.data_s1_2S
	awfn.data_p1_2P
	awfn.data_d1_2D
```

Here pp_data specifies the name and location of the file containing the 
CASINO pseudopotential. The utility then expects an input card after 
&inputpp consisting of the number of awfn.data files supplied (in this 
case 3) and then their names. The files are searched sequentially so the 
first s wavefunction found will be used for the s projector, first p for 
the p projector and so on.


*A note on the radial grid*

The utility currently performs no interpolation and attempts to use the 
same radial grid as the original pseudopotential. It therefore assumes 
that the grid will be of the standard form used by Trail and Needs.

If this is not the case the flag tn_grid=.false. can be set in the input 
file. The standard logarithmic form, r(i)=exp(xmin + i*dx) / Z is then 
assumed. Values for xmin and dx can also be specified in the input file in 
the usual way.

If interpolation from a different non-standard grid is required then the 
current recommended route is to use the casino2gon utility supplied with 
the CASINO distribution. This produces the older GON format that is 
(currently) still read by PWscf.


*Ghost states*

The Kleinman-Bylander form can unfortunately introduce ghost states into 
some calculations. If this does occur we recommend that the 
pseudopotential is re-converted using a different local channel. The local 
channel can be specified in the original CASINO pp.data file and is read 
in automatically by casino2upf.x .

### up2casino.x

This utility takes a standard UPF pseudopotential from standard input and 
writes a CASINO tabulated pseudopotential file to standard output. Any 
errors are communicated via stderr.

Usage:
	
	./up2casino.x < pseudo.UPF > pp.data

Care must be taken that the resulting pseudopotential file spec fies the 
required local channel. Also this utility should only be used with 
norm-conserving pseudopotentials.


