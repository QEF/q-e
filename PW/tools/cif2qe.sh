#!/bin/bash

#
# CIF to Quantum Espresso format converter
#               Date: 31-May-2016  Typo fixed; input file can be specified
#                                  with .cif extension (Nick Thompson)
#  Version 1.2  Date: 20-Mar-2015  Bug Fix
#  Version 1.1  Date: 23-Dec-2014  Bug Fix
#               Date: 30-Set-2014  -s option added (F. Zadra)
#  Version 1.0  Date: 15-Mar-2014  First Full conversion
#  Version 0.5  Date: 02-Oct-2013
#  Version 0.4  Date: 12 Jun 2013
#  Version 0.3  Date: 15 Nov 2012
#
# Copyright (C) 2012-2015 Carlo Nervi
# This file is distributed under the terms of the
# GNU General Public License. See the file `License'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#
# Use dos2unix to strip carriage returns at the end of the .cif files!!!!
# Tested with GNU awk v.4 - doesn't work with earlier gawk versions
#
# symmetry x,-y+1/2,z+y without ''
#

version="1.2"
USAGE="cif2qe.sh Version ${version}\nUsage: cif2qe.sh [-i] [-s] File\n    ( -i uses the ibrav of QE. -s uses CRYSTAL sg atom position mode.\n"
if [ $# == 0 -o $# -gt 2 ]; then
 printf "$USAGE"
 exit
fi

do_ibrav=0
sg=0
if [ $# == 2 ]; then
 if [ $1 == "-i" ]; then
  do_ibrav=1
  shift
 elif [ $1 == "-s" ]; then
  sg=2
  shift
 else
  printf "$USAGE"
  exit
 fi
fi

file_prefix=$(basename "$1" .cif)
file="$file_prefix".cif
if [ ! -f $file ] ; then
 echo "Error. Cannot find file $file"
 exit
fi

awk -v FILE_PREFIX="$file_prefix" -v FILE="$file" -v VERSION="$version" -v do_IBRAV=$do_ibrav -v SG=$sg '

BEGIN {
 bohr = 0.52917721092
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
# nKey             = number of recognized Keywords
# KeyW[0][1..nKey] = recognized Keywords
# KeyW[1][1..nKey] = first synonym
# ...                  ....
#

 nsynon=2
 nKey=13

 KeyW[0][1] ="_cell_length_a";               KeyW[1][1]="" 
 KeyW[0][2] ="_cell_length_b";               KeyW[1][2]=""
 KeyW[0][3] ="_cell_length_c";               KeyW[1][3]=""
 KeyW[0][4] ="_cell_angle_alpha";            KeyW[1][4]=""
 KeyW[0][5] ="_cell_angle_beta";             KeyW[1][5]=""
 KeyW[0][6] ="_cell_angle_gamma";            KeyW[1][6]=""
 KeyW[0][7] ="_atom_site_type_symbol";       KeyW[1][7]=""
 KeyW[0][8] ="_atom_site_fract_x";           KeyW[1][8]=""
 KeyW[0][9] ="_atom_site_fract_y";           KeyW[1][9]=""
 KeyW[0][10]="_atom_site_fract_z";           KeyW[1][10]=""
 KeyW[0][11]="_symmetry_equiv_pos_as_xyz";   KeyW[1][11]="_space_group_symop_operation_xyz"
 KeyW[0][12]="_symmetry_cell_setting";       KeyW[1][12]="_space_group_crystal_system"
 KeyW[0][13]="_symmetry_Int_Tables_number";  KeyW[1][13]="_space_group_IT_number"
 
#
#Tolerance for recognize identical atoms generate by symmetry
#
 tol=0.0001
#
# Separation (in A) to generate K Points
#
 separation=0.07
 totatom=0

#
# International Tables
#
# 1-2 Triclinic, 3-15 Monoclinic, 16-74 Orthorhombic, 75-142 Tetragonal, 143-167 Trigonal, 168-194 Hexagonal, 195-230 Cubic

split("P1 P-1 " \
   "P2 P2(1) C2 Pm Pc Cm Cc P2/m P2(1)/m C2/m P2/c P2(1)/c C2/c " \
   "P222 P222(1) P2(1)2(1)2 P2(1)2(1)2(1) C222(1) C222 F222 I222 I2(1)2(1)2(1) Pmm2 Pmc2(1) Pcc2 Pma2 Pca2(1) Pnc2 Pmn2(1) Pba2 Pna2(1) Pnn2 " \
       "Cmm2 Cmc2(1) Ccc2 Amm2 Abm2 Ama2 Aba2 Fmm2 Fdd2 Imm2 Iba2 Ima2 Pmmm Pnnn Pccm Pban Pmma Pnna Pmna Pcca Pbam Pccn Pbcm Pnnm Pmmn Pbcn Pbca " \
       "Pnma Cmcm Cmca Cmmm Cccm Cmma Ccca Fmmm Fddd Immm Ibam Ibca Imma " \
   "P4 P4(1) P4(2) P4(3) I4 I4(1) P-4 I-4 P4/m P4(2)/m P4/n P4(2)/n I4/m I4(1)/a P422 P42(1)2 P4(1)22 P4(1)2(1)2 P4(2)22 P4(2)2(1)2 P4(3)22 " \
       "P4(3)2(1)2 I422 I4(1)22 P4mm P4bm P4(2)cm P4(2)nm P4cc P4nc P4(2)mc P4(2)bc I4mm I4cm I4(1)md I4(1)cd P-42m P-42c P-42(1)m P-42(1)c " \
       "P-4m2 P-4c2 P-4b2 P-4n2 I-4m2 I-4c2 I-42m I-42d P4/mmm P4/mcc P4/nbm P4/nnc P4/mbm P4/mnc P4/nmm P4/ncc P4(2)/mmc P4(2)/mcm P4(2)/nbc " \
       "P4(2)/nnm P4(2)/mbc P4(2)/mnm P4(2)/nmc P4(2)/ncm I4/mmm I4/mcm I4(1)/amd I4(1)/acd " \
   "P3 P3(1) P3(2) R3 P-3 R-3 P312 P321 P3(1)12 P3(1)21 P3(2)12 P3(2)21 R32 P3m1 P31m P3c1 P31c R3m R3c P-31m P-31c P-3m1 P-3c1 R-3m R-3c " \
       "P6 P6(1) P6(5) P6(2) P6(4) P6(3) P-6 P6/m P6(3)/m P622 P6(1)22 P6(5)22 P6(2)22 P6(4)22 P6(3)22 P6mm P6cc P6(3)cm P6(3)mc P-6m2 P-6c2 P-62m " \
       "P-62c P6/mmm P6/mcc P6(3)/mcm P6(3)/mmc " \
   "P23 F23 I23 P2(1)3 I2(1)3 Pm-3 Pn-3 Fm-3 Fd-3 Im-3 Pa-3 Ia-3 P432 P4(2)32 F432 F4(1)32 I432 P4(3)32 P4(1)32 I4(1)32 P-43m F4-3m I-43m P-43n " \
       "F-43c I-43d Pm-3m Pn-3n Pm-3n Pn-3m Fm-3m Fm-3c Fd-3m Fd-3c Im-3m Ia-3d", Int_Tables, " ")
}

function parseX(str, field) {
#
# Consideriamo il testo fra due apostrofi singoli come una stringa
# ritorna il campo field-esimo della stringa str
# Se non ci sono i due apostrofi \047 ritorna il campo field-esimo
#
 idx1=0         
 str=str " "
 while (str) {
   match(str,/^ *[^ ]* |^ *\47[^\47]*\47 /)
   f=substr(str,RSTART,RLENGTH)             # save what matched in f
   gsub(/^ *|^ *\47?|\47? *$/,"",f)               # remove extra stuff
   if ( ++idx1 == field) { if (f!="") return f
     else return $field
   }
 str=substr(str,RLENGTH+1)                 # "consume" what matched
 }
 return ""
}

function search_K(str) {
 i0=0
 for (i1=1; i1<=nKey; i1++) {
   for (i2=0; i2<nsynon; i2++) {
     if (KeyW[i2][i1] == str) return i1
   }
 }
 return 0
}

function eval(cmd,expression) {
 expression = "awk \047BEGIN{printf \"%19.14f\", " cmd " }\047"
 expression |& getline l
 close(expression)
 return l     
}

#
# In str sostituisce tutte le occorrenze di pat con il valore repl es:
#   pat="x"  repl="0.945" str="x+1/2"
#
function norma(x, y, z, str) {
 gsub("x", x, str)
 gsub("y", y, str)
 gsub("z", z, str)
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

function Find_Lattice(a,b,c,alpha,beta,gamma) {
#
# find the bravais lattice name from lattice parameters
#

 thr = 1.e-4
 reticolo=""
 if ( (abs(alpha-90.0)<thr) && (abs(gamma-90.0)<thr) ) {
   if (abs(beta-90.0)<thr) {
     if ( (abs(a-b)<thr) && (abs(a-c)<thr) ) reticolo = "cubic";
     else if ( abs(a-b)<thr ) reticolo = "tetragonal";
     else reticolo = "orthorhombic";
   } else reticolo = "monoclinic"
 } else if ( (abs(alpha-90.0)<thr) && (abs(beta-90.0)<thr) && (abs(gamma-120.0)<thr) ) reticolo = "hexagonal";
   else if ( (abs(alpha-beta)<thr) && (abs(alpha-gamma)<thr) && (abs(a-b)<thr) && (abs(a-c)<thr) ) reticolo = "rhombohedral";
 else  reticolo = "triclinic"
 return reticolo
}

function Find_ibrav(spacegroup, reticolo) {
 ibrav=0
 primitive=match(spacegroup, "P")
 bodycentered=match(spacegroup, "I")
 facecentered=match(spacegroup, "F")
 basecentered=match(spacegroup, "C")
 switch (reticolo) {
 case "cubic":
   if (primitive) ibrav=1
   if (facecentered) ibrav=2
   if (bodycentered) ibrav=3
   break
 case "tetragonal":
   if (primitive) ibrav=6
   if (bodycentered) ibrav=7
   break
 case "orthorhombic":
   if (primitive) ibrav=8
   if (basecentered) ibrav=9
   if (facecentered) ibrav=10
   if (bodycentered) ibrav=11
   break;
 case "monoclinic":
   if (primitive) ibrav=-12
   if (basecentered) ibrav=13
   break
 case "triclinic":
   ibrav=14
   break
 case "hexagonal":
   ibrav=4
   break
 case "rhombohedral":
   if (primitive) ibrav=4; else ibrav=5
   break
 default: ibrav=0
 }
 return ibrav
}

function test_Var(test_type) {
 is_present=0
 for (i=0; i<Num_Var; i++) if (Var[i] == test_type) is_present=1
 return is_present
}

# Salta il commento
$1 ~ /^\#/||/^\;/ { next }

$1 ~ /loop_/ { loop_switch=1; next }

$1 ~ /^_/ {
 jvar=0
 tmp=search_K($1)
 if (loop_switch==1) {
  if (tmp > 0) Var[ivar++]=KeyW[0][tmp]; else Var[ivar++]=$1
  Num_Var=ivar
 } else {
   ivar=0
   loop_switch=0
   if (tmp > 0) Var2[KeyW[0][tmp]]=$2
 }
}

$1 ~ /^[^_]/ {
 ivar=0
 if (loop_switch==1 || loop_switch==2) {
  loop_switch=2
  for (i=0; i<Num_Var; i++) {
    Array[Var[i]][jvar]=parseX($0, i+1)
  }
  jvar++
  i = test_Var("_atom_site_type_symbol") + test_Var("_atom_site_fract_x")
  if (i == 2) natom++
  nsymm = nsymm + test_Var("_symmetry_equiv_pos_as_xyz")
 } else {
  jvar=0
  loop_switch=0
 }
}

END {
 a=strtonum(Var2["_cell_length_a"])
 b=strtonum(Var2["_cell_length_b"])
 c=strtonum(Var2["_cell_length_c"])
 alpha=strtonum(Var2["_cell_angle_alpha"])
 beta=strtonum(Var2["_cell_angle_beta"])
 gamma=strtonum(Var2["_cell_angle_gamma"])
 alphar=alpha/180.0*3.14159265358979323846
 betar=beta/180.0*3.14159265358979323846
 gammar=gamma/180.0*3.14159265358979323846
 reticolo=Find_Lattice(a,b,c,alpha,beta,gamma)
 tmpspacegroup=Var2["_symmetry_space_group"]
 tmptablenumber=Var2["_symmetry_Int_Tables_number"]
 if (tmptablenumber > 0 ) tmpspacegroup=Int_Tables[tmptablenumber]
 ibrav=Find_ibrav(tmpspacegroup, reticolo)
 KP_x = int(1./(a*separation)+0.5)
 if (KP_x < 1) KP_x=1
 KP_y = int(1./(b*separation)+0.5)
 if (KP_y < 1) KP_y=1
 KP_z = int(1./(c*separation)+0.5)
 if (KP_z < 1) KP_z=1

 ntyp=0
 for (i=0; i<natom; i++) {
  if (Type_Atom[Array["_atom_site_type_symbol"][i]] != 1) {
    Type_Atom[Array["_atom_site_type_symbol"][i]]=1
    AtomTyp[ntyp]=Array["_atom_site_type_symbol"][i]
    ntyp++
  }
 }

 for (i=0; i<natom; i++) {
  f_x=strtonum(Array["_atom_site_fract_x"][i])
  f_y=strtonum(Array["_atom_site_fract_y"][i])
  f_z=strtonum(Array["_atom_site_fract_z"][i])
  for (j=0; j<nsymm; j++) {
   if (split(Array["_symmetry_equiv_pos_as_xyz"][j], Tmp, ",") != 3) {
     print "Error in _symmetry_equiv_pos_as_xyz. Number of fields != 3: [1]=" Tmp[1] " [2]=" Tmp[2] " [3]=" Tmp[3]
     print "D: " Array["_symmetry_equiv_pos_as_xyz"][j] "  " Tmp[1] "  " Tmp[2] "  " Tmp[3]
     exit
   }
   p_X=norma(f_x, f_y, f_z, Tmp[1])
   p_Y=norma(f_x, f_y, f_z, Tmp[2])
   p_Z=norma(f_x, f_y, f_z, Tmp[3])   
   ff=0
   for (ii=0; ii<totatom; ii++) if (abs(p_X - save_X[ii]) < tol && abs(p_Y - save_Y[ii]) < tol && abs(p_Z - save_Z[ii]) < tol) ff=1
    if (ff==0) {
      save_S[totatom]=Array["_atom_site_type_symbol"][i]
      save_X[totatom]=p_X
      save_Y[totatom]=p_Y
      save_Z[totatom]=p_Z
      totatom++
    }
  }
 }
 print "! Generated by using cif2qe Version " VERSION " - Date: " strftime()
 print "!   _symmetry_space_group_name_H-M = "  Var2["_symmetry_space_group"]
 print "!   _symmetry_Int_Tables_number = " tmptablenumber
 print "!   _symmetry_cell_setting = " Var2["_symmetry_cell_setting"]
 print "! a=" a "  b=" b "  c=" c "  alpha=" alpha "  beta=" beta "  gamma=" gamma
 print "! Found by cif2qe: lattice = " reticolo "    Space group = " Int_Tables[tmptablenumber] "   ibrav = " ibrav
 print "!"
 if (SG==2)  print "! Using CRYSTALL_sg input position mode"
 print "!"
 print "! Symmetry found:"
 for (j=0; j<nsymm; j++) {
   printf "! %3i %30s    ", j+1, Array["_symmetry_equiv_pos_as_xyz"][j] 
   split(Array["_symmetry_equiv_pos_as_xyz"][j], Tmp, ",")
   printf "[%s]  [%s]  [%s]\n", Tmp[1], Tmp[2], Tmp[3]
 }
 print "&CONTROL"
 print "                       title = \x027" FILE "\x027"
 print "                 calculation = \x027" "relax\x027"
 print "                restart_mode = \x027" "from_scratch\x027"
 print "                      outdir = \x027" "./1\x027"
 print "                  pseudo_dir = \x027" "../PP/atompaw\x027"
 print "                      prefix = \x027" FILE_PREFIX "\x027"
 print "                     disk_io = \x027" "none\x027"
 print "                   verbosity = \x027" "default\x027"
 print "               etot_conv_thr = 0.0001"
 print "               forc_conv_thr = 0.001"
 print "                       nstep = 400"
 print "                     tstress = .true."
 print "                     tprnfor = .true."
 print " /"
 print " &SYSTEM"

 if (SG!=2) {
 if (do_IBRAV == 1 && ibrav != 0) {
   print "                       ibrav = " ibrav
   printf "                   celldm(1) = %19.14f\n", a/bohr
   switch (ibrav) {
    case 4:
    case 6:
    case 7:
     printf "                   celldm(3) = %19.14f\n", c/a
     break
    case 5:
    case -5:
     printf "                   celldm(4) = %19.14f\n", cos(alphar)
     break
    case 14:
     printf "                   celldm(2) = %19.14f, celldm(3) = %19.14f\n", b/a, c/a
     printf "                   celldm(4) = %19.14f\n", cos(alphar)
     printf "                   celldm(5) = %19.14f\n", cos(betar)
     printf "                   celldm(6) = %19.14f\n", cos(gammar)
     break
    case -12:
     printf "                   celldm(2) = %19.14f, celldm(3) = %19.14f\n", b/a, c/a
     printf "                   celldm(5) = %19.14f\n", cos(betar)
     break
    case 13:
     printf "                   celldm(2) = %19.14f, celldm(3) = %19.14f\n", b/a, c/a
     printf "                   celldm(4) = %19.14f\n", cos(gammar)
     break
    case 8:
    case 9:
    case 10:
    case 11:
     printf "                   celldm(2) = %19.14f, celldm(3) = %19.14f\n", b/a, c/a
     break
    case 12:
     printf "                   celldm(2) = %19.14f, celldm(3) = %19.14f\n", b/a, c/a
     printf "                   celldm(4) = %19.14f\n", cos(gammar)
     break
    case 1:
    case 2:
    case 3:
   }
 } else print "                       ibrav = 0"
 print "                         nat = " totatom
 print "                        ntyp = " ntyp
 }

if (SG==2) {
   print "                 space_group = " tmptablenumber
   print "                rhombhoedral = .TRUE."
   print "                    unique_b = .FALSE."
   print "               origin_choice = 1"
   print "                         nat = " natom
   print "                        ntyp = " ntyp
}

 print "                     ecutwfc = 60"
 print "                     ecutrho = 480"
 print "                    vdw_corr = \x027" "xdm\x027"
 print "                      xdm_a1 = 1.2153"
 print "                      xdm_a2 = 2.3704"
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
 print "\nATOMIC_SPECIES"
 for (i=0; i<ntyp; i++) printf "  %3s  %14.10f  %s.pw86pbe-n-kjpaw_psl.1.0.0.UPF\n", AtomTyp[i], Atoms[AtomTyp[i]], AtomTyp[i];

if (SG!=2) {
 print "\nATOMIC_POSITIONS crystal"
 for (i=0; i<totatom; i++)
   printf "%2s   %19.14f   %19.14f   %19.14f\n", save_S[i], save_X[i], save_Y[i], save_Z[i]
}

if (SG==2) {
 print "\nATOMIC_POSITION crystal_sg"
 for (i=0; i<natom; i++)
   printf "%2s   %19.14f    %19.14f    %19.14f\n", save_S[i], save_X[i], save_Y[i], save_Z[i]
}

 print "\nK_POINTS automatic"
 print KP_x "  " KP_y "  " KP_z "   0 0 0"

 if (do_IBRAV == 0 || SG == 2) {
   cell_px[0]=a
   cell_py[0]=0.0
   cell_pz[0]=0.0
   cell_px[1]=b*cos(gammar)
   cell_py[1]=b*sin(gammar)
   cell_pz[1]=0.0
   cell_px[2]=c*cos(betar)
   cell_py[2]=c*(cos(alphar)-cos(betar)*cos(gammar))/sin(gammar)
   cell_pz[2]=c*sqrt(1.0 - cos(alphar)^2 - cos(betar)^2 - cos(gammar)^2 + 2*cos(alphar)*cos(betar)*cos(gammar))/sin(gammar)
   print "\nCELL_PARAMETERS"
   for (i=0; i<3; i++) printf "  %19.14f   %19.14f   %19.14f\n", cell_px[i]/bohr, cell_py[i]/bohr, cell_pz[i]/bohr;
 }
 print "\n\n"

} ' $file

