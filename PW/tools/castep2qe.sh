!
! Copyright (C) 2012 Carlo Nervi
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#!/bin/bash

if [ $# != 1 ]; then
 printf "Usage: castep2qe.sh File (.cell extension not required!)\n Requires File.cell\n"
 exit
fi

if [ ! -f $1.cell ]; then
 echo "Error. Cannot find file $1.cell"
 exit
fi

#
# 1 Bohr = 0.529177249 Angstrom
#

awk -v FILE="$1" '

BEGIN {bohr = 0.52917720859;}

/%BLOCK LATTICE_CART/ {
  for (i=0; i<3; i++) {
    getline;
    cell_px[i] = $1;
    cell_py[i] = $2;
    cell_pz[i] = $3;
  }
 }

/%BLOCK POSITIONS_FRAC/ {
  natom=0;
  getline;
  do {
    frac[natom++] = $0;
    getline;
  } while ($0 !~ /%ENDBLOCK POSITIONS_FRAC/)
}

/%BLOCK KPOINTS_LIST/ {
  nks=0;
  getline;
  do {
    kpout[nks++] = sprintf ("   %19.15f   %19.15f   %19.15f     %19.15f", $1, $2, $3, $4);
    getline;
  } while ($0 !~ /%ENDBLOCK KPOINTS_LIST/)
}

/%BLOCK SPECIES_MASS/ {
  nspec=0;
  getline;
  do {
    species_symb[nspec] = $1;
    species_mass[nspec] = $2;
    nspec++;
    getline;
  } while ($0 !~ /%ENDBLOCK SPECIES_MASS/)
}

/%BLOCK SPECIES_POT/ {
  npot=0;
  getline;
  do {
    pot_symb[npot] = $1;
    pot_ps[npot] = $2;
    npot++;
    getline;
  } while ($0 !~ /%ENDBLOCK SPECIES_POT/)
}

END {
print "&CONTROL"
print "                       title = \x027" FILE "\x027"
print "                 calculation = \x027" "relax\x027"
print "                restart_mode = \x027" "from_scratch\x027"
print "                      outdir = \x27" "./1\x027"
print "                  pseudo_dir = \x027" "../PP/atompaw\x027"
print "                      prefix = \x027" FILE "\x027"
print "                     disk_io = \x027" "none\x027"
print "                   verbosity = \x027" "default\x027"
print "               etot_conv_thr = 0.00001"
print "               forc_conv_thr = 0.0001"
print "                       nstep = 680"
print "                     tstress = .true."
print "                     tprnfor = .true."
print " /"
print " &SYSTEM"
print "                       ibrav = 0"
print "                         nat = " natom
print "                        ntyp = " nspec
print "                     ecutwfc = 60"
print "                     ecutrho = 600"
print "!                      london = .true."
print " /"
print " &ELECTRONS"
print "            electron_maxstep = 200"
print "                    conv_thr = 1.0D-7"
print "              diago_thr_init = 1e-4"
print "                 startingpot = \x027" "atomic\x027"
print "                 startingwfc = \x027" "atomic\x027"
print "                 mixing_mode = \x027" "plain\x027"
print "                 mixing_beta = 0.5"
print "                 mixing_ndim = 8"
print "             diagonalization = \x027" "david\x027"
print " /"
print "&IONS"
print "                ion_dynamics = \x027" "bfgs\x027"
print " /"
print "\n"
print "ATOMIC_SPECIES"
if (npot != nspec) {
  print "Error; npot [" npot "] != nspec [" nspec "]"
  exit;
}
for (i=0; i< npot; i++) printf "  %3s  %14.10f  %s\n", species_symb[i], species_mass[i], pot_ps[i];

print "\nATOMIC_POSITIONS crystal"
for (i=0; i< natom; i++) print frac[i];

print "\nK_POINTS crystal"
print nks;
for (i=0; i<nks; i++) print kpout[i];

print "\nCELL_PARAMETERS"
for (i=0; i<3; i++) printf "  %19.15f   %19.15f   %19.15f\n", cell_px[i]/bohr, cell_py[i]/bohr, cell_pz[i]/bohr;

print "\n\n\n"

}
' $1.cell
