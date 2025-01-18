MINpuT: Structure generation for MINT calculations 


Drake Niedzielski, Eli Gerber, Yanjun Liu, Eun-Ah Kim, Tomás Arias


MINpuT is an auxiliary DFT input generation script for performing mismatched interface theory (MINT) calculations as first proposed in Ab Initio Mismatched Interface Theory of Graphene on α-RuCl3: Doping and Magnetism (2020) and later generalized into MINT-Sandwich in Niedzielski et al. (2024) [in preparation]. The script supports CIF, XSF, PWscf, POSCAR and JDFTx as output file formats for constructing a member of the MINT proxy sequence, which converges to the real incommensurate material. The output files from this script can then be used as input for standard DFT calculations. 


Mismatched interface theory (MINT) is a framework for conducting fully self-consistent electronic-structure studies of lattice-mismatched atomic heterostructures, such as graphene/α-RuCl3 or (LaSe)1.14(NbSe2)2. This approach combines the cluster approach for treating periodic structures within isolated-system methods, and the supercell approach for treating aperiodic structures within periodic methods. One material, designated as the flake, is treated as a cluster isolated from its periodic images by vacuum to prevent the introduction of unphysical strain, while the other, designated as the substrate (or slab), is kept periodic. A large enough supercell of the substrate material must be chosen to accommodate the separation of the flake. While such a flake+substrate system is only a proxy for the real incommensurate material, MINT extrapolates to the real material by studying a sequence of these proxy systems as the flake gets larger. By modeling the convergence behavior over this sequence of calculations, one can extract the target properties of the true material.  


General Notes for Running MINpuT:


MINpuT comes as either a python notebook (MINpuT.ipynb) or a python script (MINpuT.py), both require the user to make a plain-text input file (described in the next sections) which MINpuT will then parse. While the notebook creates plots on the MINT proxy system to assist with debugging, the user must manually specify the <filename>.in file by changing the respective variable at the top of the second cell in the notebook. The python script doesn’t create the same debugging plots, but the user specifies the <filename>.in file via the command line e.g.


> python MINpuT.py graphene_hBN.in 

(use "python3" on some machines)

Required packages for MINpuT.ipynb: numpy, matplotlib, re
Required packages for MINpuT.py: numpy, sys, re


CIF files for materials can easily be downloaded from the Materials Project, or similar databases. For now, MINTput only supports non-symmetrized CIF files as inputs. It is recommended that the user edit these files to contain the number of layers they desire for their flake and substrate layers as Materials project often only contains bulk layered structures with multiple layers in the unit cell.




To use the script, the user must create a plain-text input file ( <file name>.in ) containing the following parameters:


Input and Output Format:


INPUT_FORMAT: (CIF, JDFTx)


Specifies the format of the flake and substrate structure files. Currently MINpuT supports non-symmetrized CIF files and JDFTx files ( <file name>.ionpos and <file name>.lattice ).


OUTPUT_FORMAT: (CIF, JDFTx, XSF, PWscf, POSCAR)


Specifies the format in which the script will output the MINT proxy system. MINpuT’s output is intended to be for structure specification in the user’s DFT calculations. One should choose this format based on the DFT package, but CIF is the most general as there are many third party tools for converting CIF files into the input for various popular DFT packages.


Substrate and Flake Material Specification:


The user should employ the following commands to specify the materials to be used for the substrate and flake.


If INPUT_FORMAT is CIF:


SubstrateCIF: (<file name>.cif)


FlakeCIF: (<file name>.cif)


If INPUT_FORMAT is JDFTx:


SubstrateLattice: (<file name>.lattice)


SubstrateIonpos: (<file name>.ionpos)


FlakeLattice: (<file name>.lattice)


FlakeIonpos: (<file name>.ionpos)


Flake Geometry


FlakeSupercell: (nx, ny)


Used to create MINT flakes of various sizes for the convergence sequence. Specify the size of the flake’s supercell, fractional nx and ny are allowed. Further commands can be used to edit the flake geometry from this basic supercell. 


FlakeCut: (xc, yc, a, b, latticeFlake/bohr/Angstrom)


Used for trimming dangling atoms on a flake (e.g. trimming extra C on a graphene supercell so that the flake is made of hexagonal rings). Specify the center and the direction normal to a flake-cut line. The cut is performed such that flake atoms beyond (xi, yi) are removed if a(xi-xc)+b(yi-yc)<0 i.e. atoms on the positive side of the surface normal (a,b) are kept. The last argument allows the user to specify the cut in either the lattice coordinates of the flake material or cartesian coordinates with units of bohr radii or Angstrom. Multiple FlakeCut commands are allowed.


FlakeRotate: (theta, rad/deg, xc, yc, latticeFlake/bohr/Angstrom)


Used to specify flake rotation relative to the flake input file. Here, theta is the rotation angle in rad/deg. xc and yc specify the rotation center in either lattice coordinates, bohr radii, or Angstroms.


FlakeShift: (xs, ys, latticeSubstrate/latticeFlake/bohr/Angstrom)


Used to rigidly translate the flake for averaging quantities on interest over multiple interlayer registries. Here, xs and ys specify the translation vector in lattice (of either the substrate or the flake) or cartesian coordinates (in either bohr radii or Angstroms).




Flake Termination




FlakeTermination: (Y/N)


Whether to terminate the flake with additional atoms. The default terminationing species is hydrogen. For now, the user needs to edit the output files manually if they want to change the terminating species from hydrogen.


TermAtomDist: (distance, bohr/Angstrom)


Specify the distance between the terminating atoms and the atoms of the flake edge. For each edge atom, the termination(s) will lie along the directions to its truncated nearest neighbors. Terminating atoms with coordination number greater than or equal to 2 will sit directly at the truncated nearest neighbor positions and thus not be affected by this command. 




Substrate and Further Geometry Specification




MinVacuumPad (d1, d2, bohr/Angstrom)


Specify the minimum distances between atoms of adjacent flakes in the substrate supercell in either bohr radii or Angstroms. d1 and d2 refer to the distances along the two planar lattice vectors of the substrate’s unit cell lattice. The code will automatically choose the smallest supercell of the substrate that satisfies this minimum vacuum pad requirement. Using d = 0.0 can be used to enforce a commensurate lock-in in the respective direction.


InterlayerDistance (dI , bohr/Angstrom)


Specify the distance between the top substrate atoms and the bottom flake atoms in either bohr radii or Angstroms. 


LatticeVectorC (cx,cy,cz, bohr/Angstrom)


Specify the third non-planar lattice vector for the full MINT proxy system. The vector should be large enough in the z-direction to allow InterlayerDistance. If the user desires to model an isolated interface of the flake and substrate materials as in Gerber et al. (2020), cz should be chosen large enough to effectively isolate the interface its periodic images in the cz direction, as with any standard 2D ab initio calculation. By carefully choosing a smaller cz, one can model a layered 3D material in the MINT-sandwich formalism as in Niedzielski et. al. (2024) [in preparation]; cx and cy can then be used to specify the stacking registry in the full 3D material. 






Example .in file for generating MINT input configurations for twisted graphene/hBN:


________________


# MINT input (MINpuT) file that outputs structural information for a graphene flake twisted with respect to an hBN substrate


# Specify the format of the flake and substrate input files based on the DFT code
INPUT_FORMAT  CIF # JDFTx, CIF


# Specify the format that the lattice and ionpos of the MINT system should be dumped in
OUTPUT_FORMAT CIF # JDFTx, XSF, CIF, PWscf, POSCAR


# Include files to specify a unitcell of the substrate material
SubstrateCIF hBN.cif


# Include files to specify a unitcell of the flake material
FlakeCIF graphene.cif


# Specify whether to use the hydrogen termination scheme for the flake
FlakeTermination Y   #
TermAtomDist 1.09 Angstrom


# JDFTx specific hyperplane constraint commands for the flake
FlakeEdgeHPC N         # To be implemented:  Y/N add hyperplane constraints to the flake edges
FlakeXCOMHPC N         # To be implemented:  Y/N constrain the flake's X-COM
FlakeYCOMHPC N         # To be implemented:  Y/N constrain the flake's Y-COM
FlakeRotaHPC N         # To be implemented:  Y/N constrain the flake's rotation


# Specify the Flake's: size, shape, shift, and rotation  
FlakeSupercell 4 4          # can be fractional
# FlakeCut          [center-x] [center-y] [normal-x] [normal-y] [latticeFlake, bohr, Angstrom]
FlakeCut 8 0 -1 0 bohr
FlakeCut 0 2 0 1 bohr
FlakeCut 0 14 0 -1 bohr
FlakeCut 5 11.5 -1 -1 bohr
FlakeCut -1 4 1 1 bohr
FlakeCut -5 0 1 0 bohr


FlakeRotate         45 deg 1.333333 1.333333 latticeFlake   # [angle] [rad/deg] [x] [y] [latFlake, bohr, Angstrom]
FlakeShift          1.0 0.0 latticeSubstrate # [x] [y] [latFlake, latSubstrate, bohr, Angstrom]


# Used to generate a substrate supercell that accommodates the flake
MinVacuumPad 5.0 8.0 bohr # Minimum distance between periodic flake images
                            # The two numbers correspond to the substrate lattice directions
                            # 0.0 can be used to indicate a commensurate lock-in


InterlayerDistance 5.0 bohr # Distance between top of substrate and bottom of flake
LatticeVectorC 0.0 0.0 25.0 bohr # The third non-planar lattice vector
                                    # The user should think carefully about what this should be
                                    # Choose it large enough in the z-dir for a 2D system
                                    # The first two values can be used for stacking registry
________________






Upcoming features:


* Script that calls MINpuT.py that loops over flakeSupercell size to create a full sequence
* Flake hyperlane constraints in JDFTx for C.O.M. and uncooperative flake edges
