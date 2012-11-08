#!/bin/bash

#
# CIF to Quantum Espresso format converter
#  Version 0.2
#  Date: 06 Nov 2012
#
# Copyright (C) 2012 Carlo Nervi
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#

if [ $# != 1 ]; then
 printf "Usage: cif2qe.sh File (.cif extension not required!)\n Requires File.cif\n"
 exit
fi

if [ ! -f $1.cif ]; then
 echo "Error. Cannot find file $1.cif"
 exit
fi

awk -v FILE="$1" '

BEGIN {
 bohr = 0.52917720859
 nfield=split("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr " \
              "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb " \
              "Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf " \
              "Db Sg Bh Hs Mt", AtomSymb, " ")
 split("1.0079 4.0026 6.941 9.0122 10.811 12.0107 14.0067 15.9994 18.9984 20.1797 22.9897 24.305 26.9815 28.0855 30.9738 32.065 35.453 " \
       "39.948 39.0983 40.078 44.9559 47.867 50.9415 51.9961 54.938 55.845 58.9332 58.6934 63.546 65.39 69.723 72.64 74.9216 78.96 79.904 " \
       "83.8 85.4678 87.62 88.9059 91.224 92.9064 95.94 98 101.07 102.906 106.42 107.868 112.411 114.818 118.71 121.76 127.6 126.904 131.293 " \
       "132.905 137.327 138.905 140.116 140.908 144.24 145 150.36 151.964 157.25 158.925 162.5 164.93 167.259 168.934 173.04 174.967 178.49 " \
       "180.948 183.84 186.207 190.23 192.217 195.078 196.966 200.59 204.383 207.2 208.98 209 210 222 223 226 227 232.038 231.036 238.029 " \
       "237 244 243 247 247 251 252 257 258 259 262 261 262 266 264 277 268", AtomMass, " ")
 for (i=1; i<=nfield; i++) Atoms[AtomSymb[i]] = AtomMass[i]
#
#Tolerance for recognize identical atoms generate by symmetry
#
 tol=0.0001
#
# Separation (in A) to generate K Points
#
 separation=0.04
 totatom=0
}

function eval(cmd,expression) {
 expression = "awk \047BEGIN{printf \"%19.15f\", " cmd " }\047"
 expression |& getline l
 close(expression)
 return l     
}

function norma(pat, repl, str) {
 gsub(pat, repl, str)
 gsub(/--/, "+", str)
 val=eval(str)
 while (val < 0.) val+=1.
 while (val >= 1.) val-=1.
 return val
}

function abs(numero) {
 if (numero < 0) numero=-numero;
 return numero
}

# Salta il commento
/^\#/ { next
}

/loop_/ { loop_switch=1; next }


/^_/ {
 jvar=0
 if (loop_switch==1) {
  Var[ivar]=$1
  ivar++
  Num_Var=ivar
  if ($1 ~ /^_symmetry/) nsymm=0
  if ($1 ~ /^_atom/) natom=0
 } else {
   ivar=0
   loop_switch=0
   Var2[$1]=$2
 }
}

/^[^_]/ {
 ivar=0
 if (loop_switch==1 || loop_switch==2) {
  loop_switch=2
  for (i=0; i<Num_Var; i++) {
    Array[Var[i]][jvar]=$(i+1)
  }
  jvar++
  if (Var[0] ~ /^_atom/) natom++
  if (Var[0] ~ /^_symmetry/) nsymm++
 } else {
  jvar=0
  loop_switch=0
 }
}

END {
 a=strtonum(Var2["_cell_length_a"])
 b=strtonum(Var2["_cell_length_b"])
 c=strtonum(Var2["_cell_length_c"])
 alpha=strtonum(Var2["_cell_angle_alpha"])/180.0*3.14159265358979323846
 beta=strtonum(Var2["_cell_angle_beta"])/180.0*3.14159265358979323846
 gamma=strtonum(Var2["_cell_angle_gamma"])/180.0*3.14159265358979323846
 KP_x = int(1./(a*separation)+0.5)
 KP_y = int(1./(b*separation)+0.5)
 KP_z = int(1./(c*separation)+0.5)

 ntyp=0
 for (i=0; i<natom; i++) {
  if (Type_Atom[Array["_atom_site_type_symbol"][i]] != 1) {
    Type_Atom[Array["_atom_site_type_symbol"][i]]=1
    AtomTyp[ntyp]=Array["_atom_site_type_symbol"][i]
    ntyp++
  }
 }
 for (j=0; j<nsymm; j++) {
  if (split(Array["_symmetry_equiv_pos_as_xyz"][j], Tmp, ",") != 3) {
    print "Error in _symmetry_equiv_pos_as_xyz. Number of fields != 3: [1]=" Tmp[1] " [2]=" Tmp[2] " [3]=" Tmp[3]
    print "D: " Array["_symmetry_equiv_pos_as_xyz"][j] "  " Tmp[1] "  " Tmp[2] "  " Tmp[3]
    exit
  }
  TS[j][1]=Tmp[1]; TS[j][2]=Tmp[2]; TS[j][3]=Tmp[3]      
  p_X=norma("x", Array["_atom_site_fract_x"][0], Tmp[1])    
  p_Y=norma("y", Array["_atom_site_fract_y"][0], Tmp[2])
  p_Z=norma("z", Array["_atom_site_fract_z"][0], Tmp[3])
  ff=0
  for (jj=0; jj<j; jj++) {
    for (i=0; i<natom; i++) {                 
    if (abs(p_X - norma("x", Array["_atom_site_fract_x"][i], TS[jj][1])) < tol && abs(p_Y - norma("y", Array["_atom_site_fract_y"][i], TS[jj][2])) < tol && \
         abs(p_Z - norma("z", Array["_atom_site_fract_z"][i], TS[jj][3])) < tol ) ff=1
  }
  }
  if (ff==1) continue
  for (i=0; i<natom; i++) {
    save_X[totatom]=norma("x", Array["_atom_site_fract_x"][i], Tmp[1])
    save_Y[totatom]=norma("y", Array["_atom_site_fract_y"][i], Tmp[2])
    save_Z[totatom]=norma("z", Array["_atom_site_fract_z"][i], Tmp[3])
  totatom++
  }
 }
 print "&CONTROL"
 print "                       title = \x027" FILE "\x027"
# print "!   _symmetry_space_group_name_H-M = " Var2["_symmetry_space_group"]
# print "!   _symmetry_Int_Tables_number = " Var2["_symmetry_Int_Tables_number"]
# print "!   _symmetry_cell_setting = " Var2["_symmetry_cell_setting"]
 print "                 calculation = \x027" "relax\x027"
 print "                restart_mode = \x027" "from_scratch\x027"
 print "                      outdir = \x027" "./1\x027"
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
 print "                         nat = " totatom
 print "                        ntyp = " ntyp
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
 for (i=0; i<ntyp; i++) printf "  %3s  %14.10f  %s.pbe-van_ak.UPF\n", AtomTyp[i], Atoms[AtomTyp[i]], AtomTyp[i];

 print "\nATOMIC_POSITIONS crystal"
 for (i=0; i<totatom; i++)
   printf "%2s   %19.15f   %19.15f   %19.15f\n", Array["_atom_site_type_symbol"][i % natom], save_X[i], save_Y[i], save_Z[i]
 print "\nK_POINTS automatic"
 print KP_x "  " KP_y "  " KP_z "   0 0 0\n"

 cell_px[0]=a
 cell_py[0]=0.0
 cell_pz[0]=0.0
 cell_px[1]=b*cos(gamma)
 cell_py[1]=b*sin(gamma)
 cell_pz[1]=0.0
 cell_px[2]=c*cos(beta)
 cell_py[2]=c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
 cell_pz[2]=c*sqrt(1.0 - cos(alpha)^2 - cos(beta)^2 - cos(gamma)^2 + 2*cos(alpha)*cos(beta)*cos(gamma))/sin(gamma)

 print "\nCELL_PARAMETERS"
# Debug:
# print "a=" a "  b=" b "  c=" c "   alfa=" alpha "   beta=" beta "   gamma=" gamma
# for (i=0; i<3; i++) printf "  %19.15f   %19.15f   %19.15f\n", cell_px[i], cell_py[i], cell_pz[i];
 for (i=0; i<3; i++) printf "  %19.15f   %19.15f   %19.15f\n", cell_px[i]/bohr, cell_py[i]/bohr, cell_pz[i]/bohr;

 print "\n\n"

} ' $1.cif
