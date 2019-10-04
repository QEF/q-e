;; qe-funcs.el
;;
;; Copyright (C) 2016 Quantum ESPRESSO group
;; This file is distributed under the terms of the
;; GNU General Public License. See the file `License'
;; in the root directory of the present distribution,
;; or http://www.gnu.org/copyleft/gpl.txt .
;;
;; Author: Anton Kokalj (tone.kokalj at ijs.si)
;;
;; Acknowledgments:
;;
;; The implementation of qe-modes package was made possible by several
;; useful and helpful resources that are gratefully acknowledged, in
;; particular: "Mode Tutorial" of Scott Andrew Borton
;; (https://www.emacswiki.org/emacs/ModeTutorial, for indentation
;; code), "Derived Mode" and "Sample Mode" pages
;; (https://www.emacswiki.org/emacs/DerivedMode,
;; https://www.emacswiki.org/emacs/SampleMode) as well as the very
;; useful resources of Xah Lee
;; (http://ergoemacs.org/emacs/elisp_syntax_coloring.html).  Last but
;; not the least Sebastijan Peljhan is acknowledged for his work on
;; `xsf-mode' that inspired the idea of writing the qe-modes.


;; This file is not part of GNU Emacs.


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.
;;
;; This lisp script is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
;;
;; Permission is granted to distribute copies of this lisp script
;; provided the copyright notice and this permission are preserved in
;; all copies.
;;
;; You should have received a copy of the GNU General Public License
;; along with this program; if not, you can either send email to this
;; program's maintainer or write to: The Free Software Foundation,
;; Inc.; 675 Massachusetts Avenue; Cambridge, MA 02139, USA.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; send bug reports to the author (tone.kokalj at ijs.si)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;; Commentary:

;; This is the `qe-funcs.el', which defines the functions for various
;; Quantum ESPRESSO (QE) major modes. For example, functions for
;; inserting the input templates, and functions for the insertion of
;; each namelist, variable, and card. These functions follow the
;; following naming convention:
;;
;; M-x prog-insert-template
;; M-x prog-NAMELIST
;; M-x prog-CARD
;; M-X prog-variable,
;;
;; where:
;;
;; * "prog" is the lowercase name of respective program without the
;;   .x suffix (i.e. it is the lowercase variant of the PROG in the
;;   respective INPUT_PROG.html filename)
;;
;; * "NAMELIST" is the uppercase name for a given Fortran namelist
;;
;; * "CARD" is the uppercase name for a given card
;;
;; * "variable" is the lowercase name for a given namelist variable
;;
;; Note that in the above commands the spelling of namelist and card
;; names are intentionally made uppercase as to differentiate them from
;; the names of QE variables which are intentionally made lowercase.

;;; Code:



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; utility functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pw-insert-template ()
  (interactive)
  (insert "&CONTROL 
   calculation = 'relax'
/

&SYSTEM
   ! ibrav:  0 = free lattice,   1 = PC,   2 = FCC,   3 = BCC,
   !         4 = hex or trigonal P,
   !         5 = trigonal R (axis c),   -5 = trigonal R (axis <111>),
   !         6 = tetragonal P,   7 = tetragonal I,
   !         8 = orthorombic P,   9 = orthorombic base-C,   -9 = as 9 (alter description),
   !         10 = orthorombic FC,   11 = orthorombic body-C,
   !         12 = monoclinic P (axis c),   -12 = monoclinic P (axis b),
   !         13 = monoclinic base-C,   14 = triclinic
   ibrav     = 0
   celldm(1) = 1
   nat       = 1
   ntyp      = 1
   ecutwfc   = 30.0
/

&ELECTRONS
   conv_thr = 1d-7
/

&IONS
/

CELL_PARAMETERS { alat | bohr | angstrom } 
   1.00   0.00   0.00
   0.00   1.00   0.00
   0.00   0.00   1.00

ATOMIC_SPECIES
   atomLabel   atomMass   atomPseudoPotential

ATOMIC_POSITIONS { alat | bohr | angstrom | crystal | crystal_sg } 
   atomLabel   0.00   0.00   0.00

K_POINTS { tpiba | automatic | crystal | gamma | tpiba_b | crystal_b | tpiba_c | crystal_c } 
   ...insert-if-not-gamma...
"))


(defun projwfc-insert-template ()
  (interactive)
  (insert "&PROJWFC
   ! ngauss:
   !     0  = Simple Gaussian (default)
   !     1  = Methfessel-Paxton of order 1
   !    -1  = Marzari-Vanderbilt \"cold smearing\"
   !    -99 = Fermi-Dirac function
   
   ngauss  = 0 ,
   degauss = 0.007 , ! in Ry
   DeltaE  = 0.02 ,  ! in eV
   
   prefix  = 'prefix' ,
   outdir  = 'dir' ,
   filpdos = 'prefix' ,
/
"))


(defun dos-insert-template ()
  (interactive)
  (insert "&DOS
   ! ngauss:
   !     0  = Simple Gaussian (default)
   !     1  = Methfessel-Paxton of order 1
   !    -1  = Marzari-Vanderbilt \"cold smearing\"
   !    -99 = Fermi-Dirac function

   ngauss  = 0 ,
   degauss = 0.007 ,
   DeltaE  = 0.02 ,

   prefix = 'prefix' ,
   outdir = 'dir' ,
   fildos = 'prefix' ,
/
"))


(defun bands-insert-template ()
  (interactive)
  (insert "&BANDS
   prefix  = 'prefix' ,
   outdir  = 'dir' ,
   filband = 'file' ,
   
   plot_2d    = .true. ,
   lsym       = .false. ,
   no_overlap = .true. ,
/
"))


(defun pp-insert-template ()
  (interactive)
  (insert "&INPUTPP
   ! plot_num:
   !    0  = electron (pseudo-)charge density
   !    1  = total potential V_bare + V_H + V_xc
   !    2  = local ionic potential V_bare
   !    3  = local density of states at E_fermi
   !    4  = local density of electronic entropy
   !    5  = STM images
   !    6  = rho(up) - rho(down)
   !    7  = |psi|^2
   !    8  = ELF
   !    9  = rho(scf) - superposition of atomic densities
   !    10 = ILDOS
   !    11 = electrostatic potential (V_bare + V_H)
   !    12 = sawtooth electric field potential (if present)
   !    13 = noncollinear magnetization.
   !    17 = PAW all-electron valence charge density
   !    18 = XC field (noncollinear case)
   !    19 = reduced density gradient 
   !    20 = rho * second-eigenvalue-electron-density-Hessian-matrix
   !    21 = PAW all-electron charge density (valence+core).      
   plot_num = 0
   outdir   = '...'
   plot_num = 0
/

&PLOT
   nfile     = 1 
   weight(1) = 1.0

   iflag         = 3 
   output_format = 5

   fileout = '...' 
/
"))


(defun ld1-insert-template ()
  (interactive)
  (insert "&INPUT
   title =''
   
   !! either zed or atom
   zed   =   ! atomic number
   !atom =   ! atomic symbol
   
   rel = 1   ! allowed values: 0, 1, 2

   config  = '',
   iswitch =    ! allowed values: 1, 2, 3, 4
   dft     = 'PBE'
/

&INPUTP
   lpaw = 
   pseudotype =    ! allowed values: 1, 2, 3
   file_pseudopw = ''
   author = '',
   
   lloc  = 
   rcloc =
   which_augfun =     ! allowed values: 'PSQ'; for PAW: 'BESSEL', 'GAUSS' , 'BG'
   rmatch_augfun_nc = ! .true. | .false.
   nlcc =             ! .true. | .false.
   new_core_ps =      ! .true. | .false.
   rcore = 
   tm =               ! .true. | .false.
 /
nwfs  
 nls(1)    nns(1)    lls(1)  	ocs(1)     ener(1)     rcut(1)     rcutus(1)  	 [ jjs(1) ]
 nls(2)    nns(2)    lls(2)  	ocs(2)     ener(2)     rcut(2)     rcutus(2)  	 [ jjs(2) ]
 . . .
 nls(nwfs) nns(nwfs) lls(nwfs)  ocs(nwfs)  ener(nwfs)  rcut(nwfs)  rcutus(nwfs)	 [ jjs(nwfs) ]
"))


(defun neb-insert-template ()
  (interactive)
  (insert "BEGIN
BEGIN_PATH_INPUT
&PATH
   string_method  = 'neb'
   nstep_path     = 100
   num_of_images  = 10 
   opt_scheme     = 'broyden'
   CI_scheme      = 'auto'
   path_thr       = 0.05
   k_max          = 0.1
   k_min          = 0.1
   first_last_opt = .false.
   minimum_image  = .false.
   use_freezing   = .true.
   use_masses     = .false.
/
END_PATH_INPUT

BEGIN_ENGINE_INPUT
&CONTROL 
/

&SYSTEM
   ibrav     = 0
   celldm(1) = 
   nat       = 
   ntyp      = 
   ecutwfc   = 
/

&ELECTRONS
   conv_thr = 1d-7
/

&IONS
/

CELL_PARAMETERS { alat | bohr | angstrom } 
   1.00   0.00   0.00
   0.00   1.00   0.00
   0.00   0.00   1.00

ATOMIC_SPECIES
   atomLabel   atomMass   atomPseudoPotential

K_POINTS { tpiba | automatic | crystal | gamma | tpiba_b | crystal_b | tpiba_c | crystal_c } 
   ...insert-if-not-gamma...

   
BEGIN_POSITIONS
FIRST_IMAGE
ATOMIC_POSITIONS { alat | bohr | angstrom | crystal | crystal_sg } 
   atomLabel   0.00   0.00   0.00   1 1 1
   ...
   
INTERMEDIATE_IMAGE
ATOMIC_POSITIONS { alat | bohr | angstrom | crystal | crystal_sg } 
   atomLabel   0.00   0.00   0.00
   ...
   
LAST_IMAGE
ATOMIC_POSITIONS { alat | bohr | angstrom | crystal | crystal_sg } 
   atomLabel   0.00   0.00   0.00
   ...   
END_POSITIONS   
END_ENGINE_INPUT
END
"))


(defun ph-insert-template ()
  (interactive)
  (insert "title-line that is reprinted on output

&INPUTPH
   !recover = .true. | .false. 

   !amass(1) = , amass(2) = ...
   
   outdir    = ''
   prefix    = ''
   fildyn    =  ''
   
   tr2_ph    = 1d-16
   alpha_mix(1) = 0.5
   
   !nat_todo =
   
   !trans  = .true. | .false.
   !epsil  = .true. | .false.
   !lrpa   = .true. | .false.
   !lnoloc = .true. | .false.
   !lraman = .true. | .false.

   !ldisp = .true. | .false.

   !nogg  = .true. | .false.
   !asr   = .true. | .false.
/
q-point(s)-specs
nat_todo-list-of-atoms 

"))


(defun dynmat-insert-template ()
  (interactive)
  (insert "&INPUT
   fildyn = ''
   filout = ''
   filmol = ''
   filxsf = ''
   ! fileig = ''

   ! q(1) = , q(2) = , q(3) =
   ! amass(1) = , amass(2) = ...
   
   !! asr options: 'no' | 'simple' | 'crystal' | 'one-dim' | 'zero-dim' | 
   ! asr = 'no'
   
   !! for asr = 'one-dim'
   ! axis =
   
   ! lperm =
   ! lplasma =   
/

"))






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pw- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pw-CELL ()
  (interactive)
  (insert "&CELL")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun pw-CONTROL ()
  (interactive)
  (insert "&CONTROL")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun pw-ELECTRONS ()
  (interactive)
  (insert "&ELECTRONS")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun pw-IONS ()
  (interactive)
  (insert "&IONS")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun pw-SYSTEM ()
  (interactive)
  (insert "&SYSTEM")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pw- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pw-a ()
  (interactive)
  (let ((value (read-string "Value of A: ")))
    (insert "A = " value))
  )


(defun pw-adaptive_thr ()
  (interactive)
  (let ((value (read-string "Value of adaptive_thr: ")))
    (insert "adaptive_thr = " value))
  )


(defun pw-angle1 ()
  (interactive)
  (let ((value (read-string "Value of angle1: ")))
    (insert "angle1 = " value))
  )


(defun pw-angle2 ()
  (interactive)
  (let ((value (read-string "Value of angle2: ")))
    (insert "angle2 = " value))
  )


(defun pw-assume_isolated ()
  (interactive)
  (let ((value (read-string "Value of assume_isolated: ")))
    (insert "assume_isolated = '" value "'"))
  (backward-char 1)
  )


(defun pw-b ()
  (interactive)
  (let ((value (read-string "Value of B: ")))
    (insert "B = " value))
  )


(defun pw-bfgs_ndim ()
  (interactive)
  (let ((value (read-string "Value of bfgs_ndim: ")))
    (insert "bfgs_ndim = " value))
  )


(defun pw-block ()
  (interactive)
  (let ((value (read-string "Value of block: ")))
    (insert "block = " value))
  )


(defun pw-block_1 ()
  (interactive)
  (let ((value (read-string "Value of block_1: ")))
    (insert "block_1 = " value))
  )


(defun pw-block_2 ()
  (interactive)
  (let ((value (read-string "Value of block_2: ")))
    (insert "block_2 = " value))
  )


(defun pw-block_height ()
  (interactive)
  (let ((value (read-string "Value of block_height: ")))
    (insert "block_height = " value))
  )


(defun pw-c ()
  (interactive)
  (let ((value (read-string "Value of C: ")))
    (insert "C = " value))
  )


(defun pw-calculation ()
  (interactive)
  (let ((value (read-string "Value of calculation: ")))
    (insert "calculation = '" value "'"))
  (backward-char 1)
  )


(defun pw-cell_dofree ()
  (interactive)
  (let ((value (read-string "Value of cell_dofree: ")))
    (insert "cell_dofree = '" value "'"))
  (backward-char 1)
  )


(defun pw-cell_dynamics ()
  (interactive)
  (let ((value (read-string "Value of cell_dynamics: ")))
    (insert "cell_dynamics = '" value "'"))
  (backward-char 1)
  )


(defun pw-cell_factor ()
  (interactive)
  (let ((value (read-string "Value of cell_factor: ")))
    (insert "cell_factor = " value))
  )


(defun pw-celldm ()
  (interactive)
  (let ((value (read-string "Value of celldm: ")))
    (insert "celldm = " value))
  )


(defun pw-constrained_magnetization ()
  (interactive)
  (let ((value (read-string "Value of constrained_magnetization: ")))
    (insert "constrained_magnetization = '" value "'"))
  (backward-char 1)
  )


(defun pw-conv_thr ()
  (interactive)
  (let ((value (read-string "Value of conv_thr: ")))
    (insert "conv_thr = " value))
  )


(defun pw-conv_thr_init ()
  (interactive)
  (let ((value (read-string "Value of conv_thr_init: ")))
    (insert "conv_thr_init = " value))
  )


(defun pw-conv_thr_multi ()
  (interactive)
  (let ((value (read-string "Value of conv_thr_multi: ")))
    (insert "conv_thr_multi = " value))
  )


(defun pw-cosab ()
  (interactive)
  (let ((value (read-string "Value of cosAB: ")))
    (insert "cosAB = " value))
  )


(defun pw-cosac ()
  (interactive)
  (let ((value (read-string "Value of cosAC: ")))
    (insert "cosAC = " value))
  )


(defun pw-cosbc ()
  (interactive)
  (let ((value (read-string "Value of cosBC: ")))
    (insert "cosBC = " value))
  )


(defun pw-degauss ()
  (interactive)
  (let ((value (read-string "Value of degauss: ")))
    (insert "degauss = " value))
  )


(defun pw-delta_t ()
  (interactive)
  (let ((value (read-string "Value of delta_t: ")))
    (insert "delta_t = " value))
  )


(defun pw-dftd3_threebody ()
  (interactive)
  (let ((value (read-string "Value of dftd3_threebody: ")))
    (insert "dftd3_threebody = " value))
  )


(defun pw-dftd3_version ()
  (interactive)
  (let ((value (read-string "Value of dftd3_version: ")))
    (insert "dftd3_version = " value))
  )


(defun pw-diago_cg_maxiter ()
  (interactive)
  (let ((value (read-string "Value of diago_cg_maxiter: ")))
    (insert "diago_cg_maxiter = " value))
  )


(defun pw-diago_david_ndim ()
  (interactive)
  (let ((value (read-string "Value of diago_david_ndim: ")))
    (insert "diago_david_ndim = " value))
  )


(defun pw-diago_full_acc ()
  (interactive)
  (let ((value (read-string "Value of diago_full_acc: ")))
    (insert "diago_full_acc = " value))
  )


(defun pw-diago_thr_init ()
  (interactive)
  (let ((value (read-string "Value of diago_thr_init: ")))
    (insert "diago_thr_init = " value))
  )


(defun pw-diagonalization ()
  (interactive)
  (let ((value (read-string "Value of diagonalization: ")))
    (insert "diagonalization = '" value "'"))
  (backward-char 1)
  )


(defun pw-dipfield ()
  (interactive)
  (let ((value (read-string "Value of dipfield: ")))
    (insert "dipfield = " value))
  )


(defun pw-disk_io ()
  (interactive)
  (let ((value (read-string "Value of disk_io: ")))
    (insert "disk_io = '" value "'"))
  (backward-char 1)
  )


(defun pw-dt ()
  (interactive)
  (let ((value (read-string "Value of dt: ")))
    (insert "dt = " value))
  )


(defun pw-eamp ()
  (interactive)
  (let ((value (read-string "Value of eamp: ")))
    (insert "eamp = " value))
  )


(defun pw-ecfixed ()
  (interactive)
  (let ((value (read-string "Value of ecfixed: ")))
    (insert "ecfixed = " value))
  )


(defun pw-ecutfock ()
  (interactive)
  (let ((value (read-string "Value of ecutfock: ")))
    (insert "ecutfock = " value))
  )


(defun pw-ecutrho ()
  (interactive)
  (let ((value (read-string "Value of ecutrho: ")))
    (insert "ecutrho = " value))
  )


(defun pw-ecutvcut ()
  (interactive)
  (let ((value (read-string "Value of ecutvcut: ")))
    (insert "ecutvcut = " value))
  )


(defun pw-ecutwfc ()
  (interactive)
  (let ((value (read-string "Value of ecutwfc: ")))
    (insert "ecutwfc = " value))
  )


(defun pw-edir ()
  (interactive)
  (let ((value (read-string "Value of edir: ")))
    (insert "edir = " value))
  )


(defun pw-efield ()
  (interactive)
  (let ((value (read-string "Value of efield: ")))
    (insert "efield = " value))
  )


(defun pw-efield_cart ()
  (interactive)
  (let ((value (read-string "Value of efield_cart: ")))
    (insert "efield_cart = " value))
  )


(defun pw-efield_phase ()
  (interactive)
  (let ((value (read-string "Value of efield_phase: ")))
    (insert "efield_phase = '" value "'"))
  (backward-char 1)
  )


(defun pw-electron_maxstep ()
  (interactive)
  (let ((value (read-string "Value of electron_maxstep: ")))
    (insert "electron_maxstep = " value))
  )


(defun pw-emaxpos ()
  (interactive)
  (let ((value (read-string "Value of emaxpos: ")))
    (insert "emaxpos = " value))
  )


(defun pw-eopreg ()
  (interactive)
  (let ((value (read-string "Value of eopreg: ")))
    (insert "eopreg = " value))
  )


(defun pw-esm_bc ()
  (interactive)
  (let ((value (read-string "Value of esm_bc: ")))
    (insert "esm_bc = '" value "'"))
  (backward-char 1)
  )


(defun pw-esm_efield ()
  (interactive)
  (let ((value (read-string "Value of esm_efield: ")))
    (insert "esm_efield = " value))
  )


(defun pw-esm_nfit ()
  (interactive)
  (let ((value (read-string "Value of esm_nfit: ")))
    (insert "esm_nfit = " value))
  )


(defun pw-esm_w ()
  (interactive)
  (let ((value (read-string "Value of esm_w: ")))
    (insert "esm_w = " value))
  )


(defun pw-etot_conv_thr ()
  (interactive)
  (let ((value (read-string "Value of etot_conv_thr: ")))
    (insert "etot_conv_thr = " value))
  )


(defun pw-exx_fraction ()
  (interactive)
  (let ((value (read-string "Value of exx_fraction: ")))
    (insert "exx_fraction = " value))
  )


(defun pw-exxdiv_treatment ()
  (interactive)
  (let ((value (read-string "Value of exxdiv_treatment: ")))
    (insert "exxdiv_treatment = '" value "'"))
  (backward-char 1)
  )


(defun pw-fcp_mu ()
  (interactive)
  (let ((value (read-string "Value of fcp_mu: ")))
    (insert "fcp_mu = " value))
  )


(defun pw-fixed_magnetization ()
  (interactive)
  (let ((value (read-string "Value of fixed_magnetization: ")))
    (insert "fixed_magnetization = " value))
  )


(defun pw-forc_conv_thr ()
  (interactive)
  (let ((value (read-string "Value of forc_conv_thr: ")))
    (insert "forc_conv_thr = " value))
  )


(defun pw-force_symmorphic ()
  (interactive)
  (let ((value (read-string "Value of force_symmorphic: ")))
    (insert "force_symmorphic = " value))
  )


(defun pw-gate ()
  (interactive)
  (let ((value (read-string "Value of gate: ")))
    (insert "gate = " value))
  )


(defun pw-gdir ()
  (interactive)
  (let ((value (read-string "Value of gdir: ")))
    (insert "gdir = " value))
  )


(defun pw-hubbard_alpha ()
  (interactive)
  (let ((value (read-string "Value of Hubbard_alpha: ")))
    (insert "Hubbard_alpha = " value))
  )


(defun pw-hubbard_beta ()
  (interactive)
  (let ((value (read-string "Value of Hubbard_beta: ")))
    (insert "Hubbard_beta = " value))
  )


(defun pw-hubbard_j ()
  (interactive)
  (let ((value (read-string "Value of Hubbard_J: ")))
    (insert "Hubbard_J = " value))
  )


(defun pw-hubbard_j0 ()
  (interactive)
  (let ((value (read-string "Value of Hubbard_J0: ")))
    (insert "Hubbard_J0 = " value))
  )


(defun pw-hubbard_u ()
  (interactive)
  (let ((value (read-string "Value of Hubbard_U: ")))
    (insert "Hubbard_U = " value))
  )


(defun pw-ibrav ()
  (interactive)
  (let ((value (read-string "Value of ibrav: ")))
    (insert "ibrav = " value))
  )


(defun pw-input_dft ()
  (interactive)
  (let ((value (read-string "Value of input_dft: ")))
    (insert "input_dft = '" value "'"))
  (backward-char 1)
  )


(defun pw-ion_dynamics ()
  (interactive)
  (let ((value (read-string "Value of ion_dynamics: ")))
    (insert "ion_dynamics = '" value "'"))
  (backward-char 1)
  )


(defun pw-ion_positions ()
  (interactive)
  (let ((value (read-string "Value of ion_positions: ")))
    (insert "ion_positions = '" value "'"))
  (backward-char 1)
  )


(defun pw-ion_temperature ()
  (interactive)
  (let ((value (read-string "Value of ion_temperature: ")))
    (insert "ion_temperature = '" value "'"))
  (backward-char 1)
  )


(defun pw-iprint ()
  (interactive)
  (let ((value (read-string "Value of iprint: ")))
    (insert "iprint = " value))
  )


(defun pw-lambda ()
  (interactive)
  (let ((value (read-string "Value of lambda: ")))
    (insert "lambda = " value))
  )


(defun pw-lberry ()
  (interactive)
  (let ((value (read-string "Value of lberry: ")))
    (insert "lberry = " value))
  )


(defun pw-lda_plus_u ()
  (interactive)
  (let ((value (read-string "Value of lda_plus_u: ")))
    (insert "lda_plus_u = " value))
  )


(defun pw-lda_plus_u_kind ()
  (interactive)
  (let ((value (read-string "Value of lda_plus_u_kind: ")))
    (insert "lda_plus_u_kind = " value))
  )


(defun pw-lelfield ()
  (interactive)
  (let ((value (read-string "Value of lelfield: ")))
    (insert "lelfield = " value))
  )


(defun pw-lfcpopt ()
  (interactive)
  (let ((value (read-string "Value of lfcpopt: ")))
    (insert "lfcpopt = " value))
  )


(defun pw-lforcet ()
  (interactive)
  (let ((value (read-string "Value of lforcet: ")))
    (insert "lforcet = " value))
  )


(defun pw-lkpoint_dir ()
  (interactive)
  (let ((value (read-string "Value of lkpoint_dir: ")))
    (insert "lkpoint_dir = " value))
  )


(defun pw-localization_thr ()
  (interactive)
  (let ((value (read-string "Value of localization_thr: ")))
    (insert "localization_thr = " value))
  )


(defun pw-london ()
  (interactive)
  (let ((value (read-string "Value of london: ")))
    (insert "london = " value))
  )


(defun pw-london_c6 ()
  (interactive)
  (let ((value (read-string "Value of london_c6: ")))
    (insert "london_c6 = " value))
  )


(defun pw-london_rcut ()
  (interactive)
  (let ((value (read-string "Value of london_rcut: ")))
    (insert "london_rcut = " value))
  )


(defun pw-london_rvdw ()
  (interactive)
  (let ((value (read-string "Value of london_rvdw: ")))
    (insert "london_rvdw = " value))
  )


(defun pw-london_s6 ()
  (interactive)
  (let ((value (read-string "Value of london_s6: ")))
    (insert "london_s6 = " value))
  )


(defun pw-lorbm ()
  (interactive)
  (let ((value (read-string "Value of lorbm: ")))
    (insert "lorbm = " value))
  )


(defun pw-lspinorb ()
  (interactive)
  (let ((value (read-string "Value of lspinorb: ")))
    (insert "lspinorb = " value))
  )


(defun pw-max_seconds ()
  (interactive)
  (let ((value (read-string "Value of max_seconds: ")))
    (insert "max_seconds = " value))
  )


(defun pw-mixing_beta ()
  (interactive)
  (let ((value (read-string "Value of mixing_beta: ")))
    (insert "mixing_beta = " value))
  )


(defun pw-mixing_fixed_ns ()
  (interactive)
  (let ((value (read-string "Value of mixing_fixed_ns: ")))
    (insert "mixing_fixed_ns = " value))
  )


(defun pw-mixing_mode ()
  (interactive)
  (let ((value (read-string "Value of mixing_mode: ")))
    (insert "mixing_mode = '" value "'"))
  (backward-char 1)
  )


(defun pw-mixing_ndim ()
  (interactive)
  (let ((value (read-string "Value of mixing_ndim: ")))
    (insert "mixing_ndim = " value))
  )


(defun pw-nat ()
  (interactive)
  (let ((value (read-string "Value of nat: ")))
    (insert "nat = " value))
  )


(defun pw-nberrycyc ()
  (interactive)
  (let ((value (read-string "Value of nberrycyc: ")))
    (insert "nberrycyc = " value))
  )


(defun pw-nbnd ()
  (interactive)
  (let ((value (read-string "Value of nbnd: ")))
    (insert "nbnd = " value))
  )


(defun pw-no_t_rev ()
  (interactive)
  (let ((value (read-string "Value of no_t_rev: ")))
    (insert "no_t_rev = " value))
  )


(defun pw-noinv ()
  (interactive)
  (let ((value (read-string "Value of noinv: ")))
    (insert "noinv = " value))
  )


(defun pw-noncolin ()
  (interactive)
  (let ((value (read-string "Value of noncolin: ")))
    (insert "noncolin = " value))
  )


(defun pw-nosym ()
  (interactive)
  (let ((value (read-string "Value of nosym: ")))
    (insert "nosym = " value))
  )


(defun pw-nosym_evc ()
  (interactive)
  (let ((value (read-string "Value of nosym_evc: ")))
    (insert "nosym_evc = " value))
  )


(defun pw-nppstr ()
  (interactive)
  (let ((value (read-string "Value of nppstr: ")))
    (insert "nppstr = " value))
  )


(defun pw-nqx1 ()
  (interactive)
  (let ((value (read-string "Value of nqx1: ")))
    (insert "nqx1 = " value))
  )


(defun pw-nqx2 ()
  (interactive)
  (let ((value (read-string "Value of nqx2: ")))
    (insert "nqx2 = " value))
  )


(defun pw-nqx3 ()
  (interactive)
  (let ((value (read-string "Value of nqx3: ")))
    (insert "nqx3 = " value))
  )


(defun pw-nr1 ()
  (interactive)
  (let ((value (read-string "Value of nr1: ")))
    (insert "nr1 = " value))
  )


(defun pw-nr1s ()
  (interactive)
  (let ((value (read-string "Value of nr1s: ")))
    (insert "nr1s = " value))
  )


(defun pw-nr2 ()
  (interactive)
  (let ((value (read-string "Value of nr2: ")))
    (insert "nr2 = " value))
  )


(defun pw-nr2s ()
  (interactive)
  (let ((value (read-string "Value of nr2s: ")))
    (insert "nr2s = " value))
  )


(defun pw-nr3 ()
  (interactive)
  (let ((value (read-string "Value of nr3: ")))
    (insert "nr3 = " value))
  )


(defun pw-nr3s ()
  (interactive)
  (let ((value (read-string "Value of nr3s: ")))
    (insert "nr3s = " value))
  )


(defun pw-nraise ()
  (interactive)
  (let ((value (read-string "Value of nraise: ")))
    (insert "nraise = " value))
  )


(defun pw-nspin ()
  (interactive)
  (let ((value (read-string "Value of nspin: ")))
    (insert "nspin = " value))
  )


(defun pw-nstep ()
  (interactive)
  (let ((value (read-string "Value of nstep: ")))
    (insert "nstep = " value))
  )


(defun pw-ntyp ()
  (interactive)
  (let ((value (read-string "Value of ntyp: ")))
    (insert "ntyp = " value))
  )


(defun pw-occupations ()
  (interactive)
  (let ((value (read-string "Value of occupations: ")))
    (insert "occupations = '" value "'"))
  (backward-char 1)
  )


(defun pw-one_atom_occupations ()
  (interactive)
  (let ((value (read-string "Value of one_atom_occupations: ")))
    (insert "one_atom_occupations = " value))
  )


(defun pw-origin_choice ()
  (interactive)
  (let ((value (read-string "Value of origin_choice: ")))
    (insert "origin_choice = " value))
  )


(defun pw-ortho_para ()
  (interactive)
  (let ((value (read-string "Value of ortho_para: ")))
    (insert "ortho_para = " value))
  )


(defun pw-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun pw-pot_extrapolation ()
  (interactive)
  (let ((value (read-string "Value of pot_extrapolation: ")))
    (insert "pot_extrapolation = '" value "'"))
  (backward-char 1)
  )


(defun pw-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun pw-press ()
  (interactive)
  (let ((value (read-string "Value of press: ")))
    (insert "press = " value))
  )


(defun pw-press_conv_thr ()
  (interactive)
  (let ((value (read-string "Value of press_conv_thr: ")))
    (insert "press_conv_thr = " value))
  )


(defun pw-pseudo_dir ()
  (interactive)
  (let ((value (read-directory-name "Value of pseudo_dir: ")))
    (insert "pseudo_dir = '" value "'"))
  (backward-char 1)
  )


(defun pw-q2sigma ()
  (interactive)
  (let ((value (read-string "Value of q2sigma: ")))
    (insert "q2sigma = " value))
  )


(defun pw-qcutz ()
  (interactive)
  (let ((value (read-string "Value of qcutz: ")))
    (insert "qcutz = " value))
  )


(defun pw-real_space ()
  (interactive)
  (let ((value (read-string "Value of real_space: ")))
    (insert "real_space = " value))
  )


(defun pw-refold_pos ()
  (interactive)
  (let ((value (read-string "Value of refold_pos: ")))
    (insert "refold_pos = " value))
  )


(defun pw-relaxz ()
  (interactive)
  (let ((value (read-string "Value of relaxz: ")))
    (insert "relaxz = " value))
  )


(defun pw-remove_rigid_rot ()
  (interactive)
  (let ((value (read-string "Value of remove_rigid_rot: ")))
    (insert "remove_rigid_rot = " value))
  )


(defun pw-report ()
  (interactive)
  (let ((value (read-string "Value of report: ")))
    (insert "report = " value))
  )


(defun pw-restart_mode ()
  (interactive)
  (let ((value (read-string "Value of restart_mode: ")))
    (insert "restart_mode = '" value "'"))
  (backward-char 1)
  )


(defun pw-rhombohedral ()
  (interactive)
  (let ((value (read-string "Value of rhombohedral: ")))
    (insert "rhombohedral = " value))
  )


(defun pw-scf_must_converge ()
  (interactive)
  (let ((value (read-string "Value of scf_must_converge: ")))
    (insert "scf_must_converge = " value))
  )


(defun pw-screening_parameter ()
  (interactive)
  (let ((value (read-string "Value of screening_parameter: ")))
    (insert "screening_parameter = " value))
  )


(defun pw-smearing ()
  (interactive)
  (let ((value (read-string "Value of smearing: ")))
    (insert "smearing = '" value "'"))
  (backward-char 1)
  )


(defun pw-space_group ()
  (interactive)
  (let ((value (read-string "Value of space_group: ")))
    (insert "space_group = " value))
  )


(defun pw-starting_charge ()
  (interactive)
  (let ((value (read-string "Value of starting_charge: ")))
    (insert "starting_charge = " value))
  )


(defun pw-starting_magnetization ()
  (interactive)
  (let ((value (read-string "Value of starting_magnetization: ")))
    (insert "starting_magnetization = " value))
  )


(defun pw-starting_ns_eigenvalue ()
  (interactive)
  (let ((value (read-string "Value of starting_ns_eigenvalue: ")))
    (insert "starting_ns_eigenvalue = " value))
  )


(defun pw-starting_spin_angle ()
  (interactive)
  (let ((value (read-string "Value of starting_spin_angle: ")))
    (insert "starting_spin_angle = " value))
  )


(defun pw-startingpot ()
  (interactive)
  (let ((value (read-string "Value of startingpot: ")))
    (insert "startingpot = '" value "'"))
  (backward-char 1)
  )


(defun pw-startingwfc ()
  (interactive)
  (let ((value (read-string "Value of startingwfc: ")))
    (insert "startingwfc = '" value "'"))
  (backward-char 1)
  )


(defun pw-tefield ()
  (interactive)
  (let ((value (read-string "Value of tefield: ")))
    (insert "tefield = " value))
  )


(defun pw-tempw ()
  (interactive)
  (let ((value (read-string "Value of tempw: ")))
    (insert "tempw = " value))
  )


(defun pw-title ()
  (interactive)
  (let ((value (read-string "Value of title: ")))
    (insert "title = '" value "'"))
  (backward-char 1)
  )


(defun pw-tolp ()
  (interactive)
  (let ((value (read-string "Value of tolp: ")))
    (insert "tolp = " value))
  )


(defun pw-tot_charge ()
  (interactive)
  (let ((value (read-string "Value of tot_charge: ")))
    (insert "tot_charge = " value))
  )


(defun pw-tot_magnetization ()
  (interactive)
  (let ((value (read-string "Value of tot_magnetization: ")))
    (insert "tot_magnetization = " value))
  )


(defun pw-tprnfor ()
  (interactive)
  (let ((value (read-string "Value of tprnfor: ")))
    (insert "tprnfor = " value))
  )


(defun pw-tqr ()
  (interactive)
  (let ((value (read-string "Value of tqr: ")))
    (insert "tqr = " value))
  )


(defun pw-trust_radius_ini ()
  (interactive)
  (let ((value (read-string "Value of trust_radius_ini: ")))
    (insert "trust_radius_ini = " value))
  )


(defun pw-trust_radius_max ()
  (interactive)
  (let ((value (read-string "Value of trust_radius_max: ")))
    (insert "trust_radius_max = " value))
  )


(defun pw-trust_radius_min ()
  (interactive)
  (let ((value (read-string "Value of trust_radius_min: ")))
    (insert "trust_radius_min = " value))
  )


(defun pw-ts_vdw_econv_thr ()
  (interactive)
  (let ((value (read-string "Value of ts_vdw_econv_thr: ")))
    (insert "ts_vdw_econv_thr = " value))
  )


(defun pw-ts_vdw_isolated ()
  (interactive)
  (let ((value (read-string "Value of ts_vdw_isolated: ")))
    (insert "ts_vdw_isolated = " value))
  )


(defun pw-tstress ()
  (interactive)
  (let ((value (read-string "Value of tstress: ")))
    (insert "tstress = " value))
  )


(defun pw-u_projection_type ()
  (interactive)
  (let ((value (read-string "Value of U_projection_type: ")))
    (insert "U_projection_type = '" value "'"))
  (backward-char 1)
  )


(defun pw-uniqueb ()
  (interactive)
  (let ((value (read-string "Value of uniqueb: ")))
    (insert "uniqueb = " value))
  )


(defun pw-upscale ()
  (interactive)
  (let ((value (read-string "Value of upscale: ")))
    (insert "upscale = " value))
  )


(defun pw-use_all_frac ()
  (interactive)
  (let ((value (read-string "Value of use_all_frac: ")))
    (insert "use_all_frac = " value))
  )


(defun pw-vdw_corr ()
  (interactive)
  (let ((value (read-string "Value of vdw_corr: ")))
    (insert "vdw_corr = '" value "'"))
  (backward-char 1)
  )


(defun pw-verbosity ()
  (interactive)
  (let ((value (read-string "Value of verbosity: ")))
    (insert "verbosity = '" value "'"))
  (backward-char 1)
  )


(defun pw-w_1 ()
  (interactive)
  (let ((value (read-string "Value of w_1: ")))
    (insert "w_1 = " value))
  )


(defun pw-w_2 ()
  (interactive)
  (let ((value (read-string "Value of w_2: ")))
    (insert "w_2 = " value))
  )


(defun pw-wf_collect ()
  (interactive)
  (let ((value (read-string "Value of wf_collect: ")))
    (insert "wf_collect = " value))
  )


(defun pw-wfc_extrapolation ()
  (interactive)
  (let ((value (read-string "Value of wfc_extrapolation: ")))
    (insert "wfc_extrapolation = '" value "'"))
  (backward-char 1)
  )


(defun pw-wfcdir ()
  (interactive)
  (let ((value (read-directory-name "Value of wfcdir: ")))
    (insert "wfcdir = '" value "'"))
  (backward-char 1)
  )


(defun pw-wmass ()
  (interactive)
  (let ((value (read-string "Value of wmass: ")))
    (insert "wmass = " value))
  )


(defun pw-x_gamma_extrapolation ()
  (interactive)
  (let ((value (read-string "Value of x_gamma_extrapolation: ")))
    (insert "x_gamma_extrapolation = " value))
  )


(defun pw-xdm ()
  (interactive)
  (let ((value (read-string "Value of xdm: ")))
    (insert "xdm = " value))
  )


(defun pw-xdm_a1 ()
  (interactive)
  (let ((value (read-string "Value of xdm_a1: ")))
    (insert "xdm_a1 = " value))
  )


(defun pw-xdm_a2 ()
  (interactive)
  (let ((value (read-string "Value of xdm_a2: ")))
    (insert "xdm_a2 = " value))
  )


(defun pw-zgate ()
  (interactive)
  (let ((value (read-string "Value of zgate: ")))
    (insert "zgate = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pw- cards functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pw-ATOMIC_FORCES ()
 (interactive)
 (insert "ATOMIC_FORCES")
 (newline 1)
 )


(defun pw-ATOMIC_POSITIONS ()
  (interactive)
  (let ((flag (read-string "Flags: { alat | bohr | angstrom | crystal | crystal_sg } ")))
    (insert "ATOMIC_POSITIONS " flag))
  (newline 1)
  )


(defun pw-ATOMIC_SPECIES ()
 (interactive)
 (insert "ATOMIC_SPECIES")
 (newline 1)
 )


(defun pw-CELL_PARAMETERS ()
  (interactive)
  (let ((flag (read-string "Flags: { alat | bohr | angstrom } ")))
    (insert "CELL_PARAMETERS " flag))
  (newline 1)
  )


(defun pw-CONSTRAINTS ()
 (interactive)
 (insert "CONSTRAINTS")
 (newline 1)
 )


(defun pw-K_POINTS ()
  (interactive)
  (let ((flag (read-string "Flags: { tpiba | automatic | crystal | gamma | tpiba_b | crystal_b | tpiba_c | crystal_c } ")))
    (insert "K_POINTS " flag))
  (newline 1)
  )


(defun pw-OCCUPATIONS ()
 (interactive)
 (insert "OCCUPATIONS")
 (newline 1)
 )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; cppp- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun cppp-INPUTPP ()
  (interactive)
  (insert "&INPUTPP")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; cppp- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun cppp-atomic_number ()
  (interactive)
  (let ((value (read-string "Value of atomic_number: ")))
    (insert "atomic_number = " value))
  )


(defun cppp-charge_density ()
  (interactive)
  (let ((value (read-string "Value of charge_density: ")))
    (insert "charge_density = '" value "'"))
  (backward-char 1)
  )


(defun cppp-fileout ()
  (interactive)
  (let ((value (read-file-name "Value of fileout: ")))
    (insert "fileout = '" value "'"))
  (backward-char 1)
  )


(defun cppp-lbinary ()
  (interactive)
  (let ((value (read-string "Value of lbinary: ")))
    (insert "lbinary = " value))
  )


(defun cppp-lcharge ()
  (interactive)
  (let ((value (read-string "Value of lcharge: ")))
    (insert "lcharge = " value))
  )


(defun cppp-ldynamics ()
  (interactive)
  (let ((value (read-string "Value of ldynamics: ")))
    (insert "ldynamics = " value))
  )


(defun cppp-lforces ()
  (interactive)
  (let ((value (read-string "Value of lforces: ")))
    (insert "lforces = " value))
  )


(defun cppp-lpdb ()
  (interactive)
  (let ((value (read-string "Value of lpdb: ")))
    (insert "lpdb = " value))
  )


(defun cppp-lrotation ()
  (interactive)
  (let ((value (read-string "Value of lrotation: ")))
    (insert "lrotation = " value))
  )


(defun cppp-ndr ()
  (interactive)
  (let ((value (read-string "Value of ndr: ")))
    (insert "ndr = " value))
  )


(defun cppp-nframes ()
  (interactive)
  (let ((value (read-string "Value of nframes: ")))
    (insert "nframes = " value))
  )


(defun cppp-np1 ()
  (interactive)
  (let ((value (read-string "Value of np1: ")))
    (insert "np1 = " value))
  )


(defun cppp-np2 ()
  (interactive)
  (let ((value (read-string "Value of np2: ")))
    (insert "np2 = " value))
  )


(defun cppp-np3 ()
  (interactive)
  (let ((value (read-string "Value of np3: ")))
    (insert "np3 = " value))
  )


(defun cppp-ns1 ()
  (interactive)
  (let ((value (read-string "Value of ns1: ")))
    (insert "ns1 = " value))
  )


(defun cppp-ns2 ()
  (interactive)
  (let ((value (read-string "Value of ns2: ")))
    (insert "ns2 = " value))
  )


(defun cppp-ns3 ()
  (interactive)
  (let ((value (read-string "Value of ns3: ")))
    (insert "ns3 = " value))
  )


(defun cppp-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun cppp-output ()
  (interactive)
  (let ((value (read-string "Value of output: ")))
    (insert "output = '" value "'"))
  (backward-char 1)
  )


(defun cppp-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun cppp-state ()
  (interactive)
  (let ((value (read-string "Value of state: ")))
    (insert "state = '" value "'"))
  (backward-char 1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; cp- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun cp-CELL ()
  (interactive)
  (insert "&CELL")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun cp-CONTROL ()
  (interactive)
  (insert "&CONTROL")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun cp-ELECTRONS ()
  (interactive)
  (insert "&ELECTRONS")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun cp-IONS ()
  (interactive)
  (insert "&IONS")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun cp-PRESS_AI ()
  (interactive)
  (insert "&PRESS_AI")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun cp-SYSTEM ()
  (interactive)
  (insert "&SYSTEM")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun cp-WANNIER ()
  (interactive)
  (insert "&WANNIER")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; cp- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun cp-a ()
  (interactive)
  (let ((value (read-string "Value of A: ")))
    (insert "A = " value))
  )


(defun cp-abisur ()
  (interactive)
  (let ((value (read-string "Value of abisur: ")))
    (insert "abisur = " value))
  )


(defun cp-abivol ()
  (interactive)
  (let ((value (read-string "Value of abivol: ")))
    (insert "abivol = " value))
  )


(defun cp-adapt ()
  (interactive)
  (let ((value (read-string "Value of adapt: ")))
    (insert "adapt = " value))
  )


(defun cp-ampre ()
  (interactive)
  (let ((value (read-string "Value of ampre: ")))
    (insert "ampre = " value))
  )


(defun cp-amprp ()
  (interactive)
  (let ((value (read-string "Value of amprp: ")))
    (insert "amprp = " value))
  )


(defun cp-assume_isolated ()
  (interactive)
  (let ((value (read-string "Value of assume_isolated: ")))
    (insert "assume_isolated = '" value "'"))
  (backward-char 1)
  )


(defun cp-b ()
  (interactive)
  (let ((value (read-string "Value of B: ")))
    (insert "B = " value))
  )


(defun cp-c ()
  (interactive)
  (let ((value (read-string "Value of C: ")))
    (insert "C = " value))
  )


(defun cp-calculation ()
  (interactive)
  (let ((value (read-string "Value of calculation: ")))
    (insert "calculation = '" value "'"))
  (backward-char 1)
  )


(defun cp-calwf ()
  (interactive)
  (let ((value (read-string "Value of calwf: ")))
    (insert "calwf = " value))
  )


(defun cp-cell_damping ()
  (interactive)
  (let ((value (read-string "Value of cell_damping: ")))
    (insert "cell_damping = " value))
  )


(defun cp-cell_dofree ()
  (interactive)
  (let ((value (read-string "Value of cell_dofree: ")))
    (insert "cell_dofree = '" value "'"))
  (backward-char 1)
  )


(defun cp-cell_dynamics ()
  (interactive)
  (let ((value (read-string "Value of cell_dynamics: ")))
    (insert "cell_dynamics = '" value "'"))
  (backward-char 1)
  )


(defun cp-cell_factor ()
  (interactive)
  (let ((value (read-string "Value of cell_factor: ")))
    (insert "cell_factor = " value))
  )


(defun cp-cell_parameters ()
  (interactive)
  (let ((value (read-string "Value of cell_parameters: ")))
    (insert "cell_parameters = '" value "'"))
  (backward-char 1)
  )


(defun cp-cell_temperature ()
  (interactive)
  (let ((value (read-string "Value of cell_temperature: ")))
    (insert "cell_temperature = '" value "'"))
  (backward-char 1)
  )


(defun cp-cell_velocities ()
  (interactive)
  (let ((value (read-string "Value of cell_velocities: ")))
    (insert "cell_velocities = '" value "'"))
  (backward-char 1)
  )


(defun cp-celldm ()
  (interactive)
  (let ((value (read-string "Value of celldm: ")))
    (insert "celldm = " value))
  )


(defun cp-conv_thr ()
  (interactive)
  (let ((value (read-string "Value of conv_thr: ")))
    (insert "conv_thr = " value))
  )


(defun cp-cosab ()
  (interactive)
  (let ((value (read-string "Value of cosAB: ")))
    (insert "cosAB = " value))
  )


(defun cp-cosac ()
  (interactive)
  (let ((value (read-string "Value of cosAC: ")))
    (insert "cosAC = " value))
  )


(defun cp-cosbc ()
  (interactive)
  (let ((value (read-string "Value of cosBC: ")))
    (insert "cosBC = " value))
  )


(defun cp-degauss ()
  (interactive)
  (let ((value (read-string "Value of degauss: ")))
    (insert "degauss = " value))
  )


(defun cp-disk_io ()
  (interactive)
  (let ((value (read-string "Value of disk_io: ")))
    (insert "disk_io = '" value "'"))
  (backward-char 1)
  )


(defun cp-dt ()
  (interactive)
  (let ((value (read-string "Value of dt: ")))
    (insert "dt = " value))
  )


(defun cp-dthr ()
  (interactive)
  (let ((value (read-string "Value of dthr: ")))
    (insert "dthr = " value))
  )


(defun cp-ecfixed ()
  (interactive)
  (let ((value (read-string "Value of ecfixed: ")))
    (insert "ecfixed = " value))
  )


(defun cp-ecutrho ()
  (interactive)
  (let ((value (read-string "Value of ecutrho: ")))
    (insert "ecutrho = " value))
  )


(defun cp-ecutwfc ()
  (interactive)
  (let ((value (read-string "Value of ecutwfc: ")))
    (insert "ecutwfc = " value))
  )


(defun cp-efield ()
  (interactive)
  (let ((value (read-string "Value of efield: ")))
    (insert "efield = " value))
  )


(defun cp-efx0 ()
  (interactive)
  (let ((value (read-string "Value of efx0: ")))
    (insert "efx0 = " value))
  )


(defun cp-efx1 ()
  (interactive)
  (let ((value (read-string "Value of efx1: ")))
    (insert "efx1 = " value))
  )


(defun cp-efy0 ()
  (interactive)
  (let ((value (read-string "Value of efy0: ")))
    (insert "efy0 = " value))
  )


(defun cp-efy1 ()
  (interactive)
  (let ((value (read-string "Value of efy1: ")))
    (insert "efy1 = " value))
  )


(defun cp-efz0 ()
  (interactive)
  (let ((value (read-string "Value of efz0: ")))
    (insert "efz0 = " value))
  )


(defun cp-efz1 ()
  (interactive)
  (let ((value (read-string "Value of efz1: ")))
    (insert "efz1 = " value))
  )


(defun cp-ekin_conv_thr ()
  (interactive)
  (let ((value (read-string "Value of ekin_conv_thr: ")))
    (insert "ekin_conv_thr = " value))
  )


(defun cp-ekincw ()
  (interactive)
  (let ((value (read-string "Value of ekincw: ")))
    (insert "ekincw = " value))
  )


(defun cp-electron_damping ()
  (interactive)
  (let ((value (read-string "Value of electron_damping: ")))
    (insert "electron_damping = " value))
  )


(defun cp-electron_dynamics ()
  (interactive)
  (let ((value (read-string "Value of electron_dynamics: ")))
    (insert "electron_dynamics = '" value "'"))
  (backward-char 1)
  )


(defun cp-electron_maxstep ()
  (interactive)
  (let ((value (read-string "Value of electron_maxstep: ")))
    (insert "electron_maxstep = " value))
  )


(defun cp-electron_temperature ()
  (interactive)
  (let ((value (read-string "Value of electron_temperature: ")))
    (insert "electron_temperature = '" value "'"))
  (backward-char 1)
  )


(defun cp-electron_velocities ()
  (interactive)
  (let ((value (read-string "Value of electron_velocities: ")))
    (insert "electron_velocities = '" value "'"))
  (backward-char 1)
  )


(defun cp-emass ()
  (interactive)
  (let ((value (read-string "Value of emass: ")))
    (insert "emass = " value))
  )


(defun cp-emass_cutoff ()
  (interactive)
  (let ((value (read-string "Value of emass_cutoff: ")))
    (insert "emass_cutoff = " value))
  )


(defun cp-epol ()
  (interactive)
  (let ((value (read-string "Value of epol: ")))
    (insert "epol = " value))
  )


(defun cp-etot_conv_thr ()
  (interactive)
  (let ((value (read-string "Value of etot_conv_thr: ")))
    (insert "etot_conv_thr = " value))
  )


(defun cp-exx_dis_cutoff ()
  (interactive)
  (let ((value (read-string "Value of exx_dis_cutoff: ")))
    (insert "exx_dis_cutoff = " value))
  )


(defun cp-exx_fraction ()
  (interactive)
  (let ((value (read-string "Value of exx_fraction: ")))
    (insert "exx_fraction = " value))
  )


(defun cp-exx_me_rcut_pair ()
  (interactive)
  (let ((value (read-string "Value of exx_me_rcut_pair: ")))
    (insert "exx_me_rcut_pair = " value))
  )


(defun cp-exx_me_rcut_self ()
  (interactive)
  (let ((value (read-string "Value of exx_me_rcut_self: ")))
    (insert "exx_me_rcut_self = " value))
  )


(defun cp-exx_neigh ()
  (interactive)
  (let ((value (read-string "Value of exx_neigh: ")))
    (insert "exx_neigh = " value))
  )


(defun cp-exx_poisson_eps ()
  (interactive)
  (let ((value (read-string "Value of exx_poisson_eps: ")))
    (insert "exx_poisson_eps = " value))
  )


(defun cp-exx_ps_rcut_pair ()
  (interactive)
  (let ((value (read-string "Value of exx_ps_rcut_pair: ")))
    (insert "exx_ps_rcut_pair = " value))
  )


(defun cp-exx_ps_rcut_self ()
  (interactive)
  (let ((value (read-string "Value of exx_ps_rcut_self: ")))
    (insert "exx_ps_rcut_self = " value))
  )


(defun cp-fnhscl ()
  (interactive)
  (let ((value (read-string "Value of fnhscl: ")))
    (insert "fnhscl = " value))
  )


(defun cp-fnosee ()
  (interactive)
  (let ((value (read-string "Value of fnosee: ")))
    (insert "fnosee = " value))
  )


(defun cp-fnoseh ()
  (interactive)
  (let ((value (read-string "Value of fnoseh: ")))
    (insert "fnoseh = " value))
  )


(defun cp-fnosep ()
  (interactive)
  (let ((value (read-string "Value of fnosep: ")))
    (insert "fnosep = " value))
  )


(defun cp-forc_conv_thr ()
  (interactive)
  (let ((value (read-string "Value of forc_conv_thr: ")))
    (insert "forc_conv_thr = " value))
  )


(defun cp-grease ()
  (interactive)
  (let ((value (read-string "Value of grease: ")))
    (insert "grease = " value))
  )


(defun cp-greash ()
  (interactive)
  (let ((value (read-string "Value of greash: ")))
    (insert "greash = " value))
  )


(defun cp-greasp ()
  (interactive)
  (let ((value (read-string "Value of greasp: ")))
    (insert "greasp = " value))
  )


(defun cp-hubbard_u ()
  (interactive)
  (let ((value (read-string "Value of Hubbard_U: ")))
    (insert "Hubbard_U = " value))
  )


(defun cp-ibrav ()
  (interactive)
  (let ((value (read-string "Value of ibrav: ")))
    (insert "ibrav = " value))
  )


(defun cp-iesr ()
  (interactive)
  (let ((value (read-string "Value of iesr: ")))
    (insert "iesr = " value))
  )


(defun cp-input_dft ()
  (interactive)
  (let ((value (read-string "Value of input_dft: ")))
    (insert "input_dft = '" value "'"))
  (backward-char 1)
  )


(defun cp-ion_damping ()
  (interactive)
  (let ((value (read-string "Value of ion_damping: ")))
    (insert "ion_damping = " value))
  )


(defun cp-ion_dynamics ()
  (interactive)
  (let ((value (read-string "Value of ion_dynamics: ")))
    (insert "ion_dynamics = '" value "'"))
  (backward-char 1)
  )


(defun cp-ion_nstepe ()
  (interactive)
  (let ((value (read-string "Value of ion_nstepe: ")))
    (insert "ion_nstepe = " value))
  )


(defun cp-ion_positions ()
  (interactive)
  (let ((value (read-string "Value of ion_positions: ")))
    (insert "ion_positions = '" value "'"))
  (backward-char 1)
  )


(defun cp-ion_radius ()
  (interactive)
  (let ((value (read-string "Value of ion_radius: ")))
    (insert "ion_radius = " value))
  )


(defun cp-ion_temperature ()
  (interactive)
  (let ((value (read-string "Value of ion_temperature: ")))
    (insert "ion_temperature = '" value "'"))
  (backward-char 1)
  )


(defun cp-ion_velocities ()
  (interactive)
  (let ((value (read-string "Value of ion_velocities: ")))
    (insert "ion_velocities = '" value "'"))
  (backward-char 1)
  )


(defun cp-iprint ()
  (interactive)
  (let ((value (read-string "Value of iprint: ")))
    (insert "iprint = " value))
  )


(defun cp-isave ()
  (interactive)
  (let ((value (read-string "Value of isave: ")))
    (insert "isave = " value))
  )


(defun cp-lambda_cold ()
  (interactive)
  (let ((value (read-string "Value of lambda_cold: ")))
    (insert "lambda_cold = " value))
  )


(defun cp-lda_plus_u ()
  (interactive)
  (let ((value (read-string "Value of lda_plus_u: ")))
    (insert "lda_plus_u = " value))
  )


(defun cp-london_rcut ()
  (interactive)
  (let ((value (read-string "Value of london_rcut: ")))
    (insert "london_rcut = " value))
  )


(defun cp-london_s6 ()
  (interactive)
  (let ((value (read-string "Value of london_s6: ")))
    (insert "london_s6 = " value))
  )


(defun cp-max_seconds ()
  (interactive)
  (let ((value (read-string "Value of max_seconds: ")))
    (insert "max_seconds = " value))
  )


(defun cp-maxiter ()
  (interactive)
  (let ((value (read-string "Value of maxiter: ")))
    (insert "maxiter = " value))
  )


(defun cp-maxwfdt ()
  (interactive)
  (let ((value (read-string "Value of maxwfdt: ")))
    (insert "maxwfdt = " value))
  )


(defun cp-memory ()
  (interactive)
  (let ((value (read-string "Value of memory: ")))
    (insert "memory = '" value "'"))
  (backward-char 1)
  )


(defun cp-n_inner ()
  (interactive)
  (let ((value (read-string "Value of n_inner: ")))
    (insert "n_inner = " value))
  )


(defun cp-nat ()
  (interactive)
  (let ((value (read-string "Value of nat: ")))
    (insert "nat = " value))
  )


(defun cp-nbnd ()
  (interactive)
  (let ((value (read-string "Value of nbnd: ")))
    (insert "nbnd = " value))
  )


(defun cp-ndega ()
  (interactive)
  (let ((value (read-string "Value of ndega: ")))
    (insert "ndega = " value))
  )


(defun cp-ndr ()
  (interactive)
  (let ((value (read-string "Value of ndr: ")))
    (insert "ndr = " value))
  )


(defun cp-ndw ()
  (interactive)
  (let ((value (read-string "Value of ndw: ")))
    (insert "ndw = " value))
  )


(defun cp-nhgrp ()
  (interactive)
  (let ((value (read-string "Value of nhgrp: ")))
    (insert "nhgrp = " value))
  )


(defun cp-nhpcl ()
  (interactive)
  (let ((value (read-string "Value of nhpcl: ")))
    (insert "nhpcl = " value))
  )


(defun cp-nhptyp ()
  (interactive)
  (let ((value (read-string "Value of nhptyp: ")))
    (insert "nhptyp = " value))
  )


(defun cp-ninter_cold_restart ()
  (interactive)
  (let ((value (read-string "Value of ninter_cold_restart: ")))
    (insert "ninter_cold_restart = " value))
  )


(defun cp-nit ()
  (interactive)
  (let ((value (read-string "Value of nit: ")))
    (insert "nit = " value))
  )


(defun cp-niter_cg_restart ()
  (interactive)
  (let ((value (read-string "Value of niter_cg_restart: ")))
    (insert "niter_cg_restart = " value))
  )


(defun cp-nr1 ()
  (interactive)
  (let ((value (read-string "Value of nr1: ")))
    (insert "nr1 = " value))
  )


(defun cp-nr1b ()
  (interactive)
  (let ((value (read-string "Value of nr1b: ")))
    (insert "nr1b = " value))
  )


(defun cp-nr1s ()
  (interactive)
  (let ((value (read-string "Value of nr1s: ")))
    (insert "nr1s = " value))
  )


(defun cp-nr2 ()
  (interactive)
  (let ((value (read-string "Value of nr2: ")))
    (insert "nr2 = " value))
  )


(defun cp-nr2b ()
  (interactive)
  (let ((value (read-string "Value of nr2b: ")))
    (insert "nr2b = " value))
  )


(defun cp-nr2s ()
  (interactive)
  (let ((value (read-string "Value of nr2s: ")))
    (insert "nr2s = " value))
  )


(defun cp-nr3 ()
  (interactive)
  (let ((value (read-string "Value of nr3: ")))
    (insert "nr3 = " value))
  )


(defun cp-nr3b ()
  (interactive)
  (let ((value (read-string "Value of nr3b: ")))
    (insert "nr3b = " value))
  )


(defun cp-nr3s ()
  (interactive)
  (let ((value (read-string "Value of nr3s: ")))
    (insert "nr3s = " value))
  )


(defun cp-nsd ()
  (interactive)
  (let ((value (read-string "Value of nsd: ")))
    (insert "nsd = " value))
  )


(defun cp-nspin ()
  (interactive)
  (let ((value (read-string "Value of nspin: ")))
    (insert "nspin = " value))
  )


(defun cp-nstep ()
  (interactive)
  (let ((value (read-string "Value of nstep: ")))
    (insert "nstep = " value))
  )


(defun cp-nsteps ()
  (interactive)
  (let ((value (read-string "Value of nsteps: ")))
    (insert "nsteps = " value))
  )


(defun cp-ntyp ()
  (interactive)
  (let ((value (read-string "Value of ntyp: ")))
    (insert "ntyp = " value))
  )


(defun cp-nwf ()
  (interactive)
  (let ((value (read-string "Value of nwf: ")))
    (insert "nwf = " value))
  )


(defun cp-occupations ()
  (interactive)
  (let ((value (read-string "Value of occupations: ")))
    (insert "occupations = '" value "'"))
  (backward-char 1)
  )


(defun cp-ortho_eps ()
  (interactive)
  (let ((value (read-string "Value of ortho_eps: ")))
    (insert "ortho_eps = " value))
  )


(defun cp-ortho_max ()
  (interactive)
  (let ((value (read-string "Value of ortho_max: ")))
    (insert "ortho_max = " value))
  )


(defun cp-ortho_para ()
  (interactive)
  (let ((value (read-string "Value of ortho_para: ")))
    (insert "ortho_para = " value))
  )


(defun cp-orthogonalization ()
  (interactive)
  (let ((value (read-string "Value of orthogonalization: ")))
    (insert "orthogonalization = '" value "'"))
  (backward-char 1)
  )


(defun cp-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun cp-p_ext ()
  (interactive)
  (let ((value (read-string "Value of P_ext: ")))
    (insert "P_ext = " value))
  )


(defun cp-p_fin ()
  (interactive)
  (let ((value (read-string "Value of P_fin: ")))
    (insert "P_fin = " value))
  )


(defun cp-p_in ()
  (interactive)
  (let ((value (read-string "Value of P_in: ")))
    (insert "P_in = " value))
  )


(defun cp-passop ()
  (interactive)
  (let ((value (read-string "Value of passop: ")))
    (insert "passop = " value))
  )


(defun cp-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun cp-press ()
  (interactive)
  (let ((value (read-string "Value of press: ")))
    (insert "press = " value))
  )


(defun cp-pseudo_dir ()
  (interactive)
  (let ((value (read-directory-name "Value of pseudo_dir: ")))
    (insert "pseudo_dir = '" value "'"))
  (backward-char 1)
  )


(defun cp-pvar ()
  (interactive)
  (let ((value (read-string "Value of pvar: ")))
    (insert "pvar = " value))
  )


(defun cp-q2sigma ()
  (interactive)
  (let ((value (read-string "Value of q2sigma: ")))
    (insert "q2sigma = " value))
  )


(defun cp-qcutz ()
  (interactive)
  (let ((value (read-string "Value of qcutz: ")))
    (insert "qcutz = " value))
  )


(defun cp-remove_rigid_rot ()
  (interactive)
  (let ((value (read-string "Value of remove_rigid_rot: ")))
    (insert "remove_rigid_rot = " value))
  )


(defun cp-restart_mode ()
  (interactive)
  (let ((value (read-string "Value of restart_mode: ")))
    (insert "restart_mode = '" value "'"))
  (backward-char 1)
  )


(defun cp-rho_thr ()
  (interactive)
  (let ((value (read-string "Value of rho_thr: ")))
    (insert "rho_thr = " value))
  )


(defun cp-saverho ()
  (interactive)
  (let ((value (read-string "Value of saverho: ")))
    (insert "saverho = " value))
  )


(defun cp-smearing ()
  (interactive)
  (let ((value (read-string "Value of smearing: ")))
    (insert "smearing = '" value "'"))
  (backward-char 1)
  )


(defun cp-startingwfc ()
  (interactive)
  (let ((value (read-string "Value of startingwfc: ")))
    (insert "startingwfc = '" value "'"))
  (backward-char 1)
  )


(defun cp-surf_t ()
  (interactive)
  (let ((value (read-string "Value of Surf_t: ")))
    (insert "Surf_t = " value))
  )


(defun cp-sw_len ()
  (interactive)
  (let ((value (read-string "Value of sw_len: ")))
    (insert "sw_len = " value))
  )


(defun cp-tabps ()
  (interactive)
  (let ((value (read-string "Value of tabps: ")))
    (insert "tabps = " value))
  )


(defun cp-tcg ()
  (interactive)
  (let ((value (read-string "Value of tcg: ")))
    (insert "tcg = " value))
  )


(defun cp-tefield ()
  (interactive)
  (let ((value (read-string "Value of tefield: ")))
    (insert "tefield = " value))
  )


(defun cp-temph ()
  (interactive)
  (let ((value (read-string "Value of temph: ")))
    (insert "temph = " value))
  )


(defun cp-tempw ()
  (interactive)
  (let ((value (read-string "Value of tempw: ")))
    (insert "tempw = " value))
  )


(defun cp-title ()
  (interactive)
  (let ((value (read-string "Value of title: ")))
    (insert "title = '" value "'"))
  (backward-char 1)
  )


(defun cp-tolp ()
  (interactive)
  (let ((value (read-string "Value of tolp: ")))
    (insert "tolp = " value))
  )


(defun cp-tolw ()
  (interactive)
  (let ((value (read-string "Value of tolw: ")))
    (insert "tolw = " value))
  )


(defun cp-tot_charge ()
  (interactive)
  (let ((value (read-string "Value of tot_charge: ")))
    (insert "tot_charge = " value))
  )


(defun cp-tot_magnetization ()
  (interactive)
  (let ((value (read-string "Value of tot_magnetization: ")))
    (insert "tot_magnetization = " value))
  )


(defun cp-tprnfor ()
  (interactive)
  (let ((value (read-string "Value of tprnfor: ")))
    (insert "tprnfor = " value))
  )


(defun cp-tranp ()
  (interactive)
  (let ((value (read-string "Value of tranp: ")))
    (insert "tranp = " value))
  )


(defun cp-ts_vdw ()
  (interactive)
  (let ((value (read-string "Value of ts_vdw: ")))
    (insert "ts_vdw = " value))
  )


(defun cp-ts_vdw_econv_thr ()
  (interactive)
  (let ((value (read-string "Value of ts_vdw_econv_thr: ")))
    (insert "ts_vdw_econv_thr = " value))
  )


(defun cp-ts_vdw_isolated ()
  (interactive)
  (let ((value (read-string "Value of ts_vdw_isolated: ")))
    (insert "ts_vdw_isolated = " value))
  )


(defun cp-tstress ()
  (interactive)
  (let ((value (read-string "Value of tstress: ")))
    (insert "tstress = " value))
  )


(defun cp-vdw_corr ()
  (interactive)
  (let ((value (read-string "Value of vdw_corr: ")))
    (insert "vdw_corr = '" value "'"))
  (backward-char 1)
  )


(defun cp-verbosity ()
  (interactive)
  (let ((value (read-string "Value of verbosity: ")))
    (insert "verbosity = '" value "'"))
  (backward-char 1)
  )


(defun cp-wf_efield ()
  (interactive)
  (let ((value (read-string "Value of wf_efield: ")))
    (insert "wf_efield = " value))
  )


(defun cp-wf_friction ()
  (interactive)
  (let ((value (read-string "Value of wf_friction: ")))
    (insert "wf_friction = " value))
  )


(defun cp-wf_q ()
  (interactive)
  (let ((value (read-string "Value of wf_q: ")))
    (insert "wf_q = " value))
  )


(defun cp-wf_switch ()
  (interactive)
  (let ((value (read-string "Value of wf_switch: ")))
    (insert "wf_switch = " value))
  )


(defun cp-wfdt ()
  (interactive)
  (let ((value (read-string "Value of wfdt: ")))
    (insert "wfdt = " value))
  )


(defun cp-wffort ()
  (interactive)
  (let ((value (read-string "Value of wffort: ")))
    (insert "wffort = " value))
  )


(defun cp-wfsd ()
  (interactive)
  (let ((value (read-string "Value of wfsd: ")))
    (insert "wfsd = " value))
  )


(defun cp-wmass ()
  (interactive)
  (let ((value (read-string "Value of wmass: ")))
    (insert "wmass = " value))
  )


(defun cp-writev ()
  (interactive)
  (let ((value (read-string "Value of writev: ")))
    (insert "writev = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; cp- cards functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun cp-ATOMIC_FORCES ()
 (interactive)
 (insert "ATOMIC_FORCES")
 (newline 1)
 )


(defun cp-ATOMIC_POSITIONS ()
  (interactive)
  (let ((flag (read-string "Flags: { alat | bohr | angstrom | crystal } ")))
    (insert "ATOMIC_POSITIONS " flag))
  (newline 1)
  )


(defun cp-ATOMIC_SPECIES ()
 (interactive)
 (insert "ATOMIC_SPECIES")
 (newline 1)
 )


(defun cp-ATOMIC_VELOCITIES ()
  (interactive)
  (let ((flag (read-string "Flags: { a.u } ")))
    (insert "ATOMIC_VELOCITIES " flag))
  (newline 1)
  )


(defun cp-AUTOPILOT ()
 (interactive)
 (insert "AUTOPILOT")
 (newline 1)
 )


(defun cp-CELL_PARAMETERS ()
  (interactive)
  (let ((flag (read-string "Flags: { bohr | angstrom | alat } ")))
    (insert "CELL_PARAMETERS " flag))
  (newline 1)
  )


(defun cp-CONSTRAINTS ()
 (interactive)
 (insert "CONSTRAINTS")
 (newline 1)
 )


(defun cp-OCCUPATIONS ()
 (interactive)
 (insert "OCCUPATIONS")
 (newline 1)
 )


(defun cp-PLOT_WANNIER ()
 (interactive)
 (insert "PLOT_WANNIER")
 (newline 1)
 )


(defun cp-REF_CELL_PARAMETERS ()
  (interactive)
  (let ((flag (read-string "Flags: { bohr | angstrom } ")))
    (insert "REF_CELL_PARAMETERS " flag))
  (newline 1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; ppacf- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ppacf-PPACF ()
  (interactive)
  (insert "&PPACF")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; ppacf- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ppacf-code_num ()
  (interactive)
  (let ((value (read-string "Value of code_num: ")))
    (insert "code_num = " value))
  )


(defun ppacf-lfock ()
  (interactive)
  (let ((value (read-string "Value of lfock: ")))
    (insert "lfock = " value))
  )


(defun ppacf-lplot ()
  (interactive)
  (let ((value (read-string "Value of lplot: ")))
    (insert "lplot = " value))
  )


(defun ppacf-ltks ()
  (interactive)
  (let ((value (read-string "Value of ltks: ")))
    (insert "ltks = " value))
  )


(defun ppacf-n_lambda ()
  (interactive)
  (let ((value (read-string "Value of n_lambda: ")))
    (insert "n_lambda = " value))
  )


(defun ppacf-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun ppacf-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun ppacf-pseudo_dir ()
  (interactive)
  (let ((value (read-directory-name "Value of pseudo_dir: ")))
    (insert "pseudo_dir = '" value "'"))
  (backward-char 1)
  )


(defun ppacf-use_ace ()
  (interactive)
  (let ((value (read-string "Value of use_ace: ")))
    (insert "use_ace = " value))
  )


(defun ppacf-vdw_table_name ()
  (interactive)
  (let ((value (read-string "Value of vdw_table_name: ")))
    (insert "vdw_table_name = '" value "'"))
  (backward-char 1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; importexport_binary- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun importexport_binary-INPUTPP ()
  (interactive)
  (insert "&INPUTPP")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; importexport_binary- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun importexport_binary-direction ()
  (interactive)
  (let ((value (read-string "Value of direction: ")))
    (insert "direction = '" value "'"))
  (backward-char 1)
  )


(defun importexport_binary-newoutdir ()
  (interactive)
  (let ((value (read-directory-name "Value of newoutdir: ")))
    (insert "newoutdir = '" value "'"))
  (backward-char 1)
  )


(defun importexport_binary-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun importexport_binary-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; projwfc- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun projwfc-PROJWFC ()
  (interactive)
  (insert "&PROJWFC")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; projwfc- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun projwfc-degauss ()
  (interactive)
  (let ((value (read-string "Value of degauss: ")))
    (insert "degauss = " value))
  )


(defun projwfc-deltae ()
  (interactive)
  (let ((value (read-string "Value of DeltaE: ")))
    (insert "DeltaE = " value))
  )


(defun projwfc-emax ()
  (interactive)
  (let ((value (read-string "Value of Emax: ")))
    (insert "Emax = " value))
  )


(defun projwfc-emin ()
  (interactive)
  (let ((value (read-string "Value of Emin: ")))
    (insert "Emin = " value))
  )


(defun projwfc-filpdos ()
  (interactive)
  (let ((value (read-string "Value of filpdos: ")))
    (insert "filpdos = '" value "'"))
  (backward-char 1)
  )


(defun projwfc-filproj ()
  (interactive)
  (let ((value (read-string "Value of filproj: ")))
    (insert "filproj = '" value "'"))
  (backward-char 1)
  )


(defun projwfc-irmax ()
  (interactive)
  (let ((value (read-string "Value of irmax: ")))
    (insert "irmax = " value))
  )


(defun projwfc-irmin ()
  (interactive)
  (let ((value (read-string "Value of irmin: ")))
    (insert "irmin = " value))
  )


(defun projwfc-kresolveddos ()
  (interactive)
  (let ((value (read-string "Value of kresolveddos: ")))
    (insert "kresolveddos = " value))
  )


(defun projwfc-lbinary_data ()
  (interactive)
  (let ((value (read-string "Value of lbinary_data: ")))
    (insert "lbinary_data = " value))
  )


(defun projwfc-lsym ()
  (interactive)
  (let ((value (read-string "Value of lsym: ")))
    (insert "lsym = " value))
  )


(defun projwfc-lwrite_overlaps ()
  (interactive)
  (let ((value (read-string "Value of lwrite_overlaps: ")))
    (insert "lwrite_overlaps = " value))
  )


(defun projwfc-n_proj_boxes ()
  (interactive)
  (let ((value (read-string "Value of n_proj_boxes: ")))
    (insert "n_proj_boxes = " value))
  )


(defun projwfc-ngauss ()
  (interactive)
  (let ((value (read-string "Value of ngauss: ")))
    (insert "ngauss = " value))
  )


(defun projwfc-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun projwfc-pawproj ()
  (interactive)
  (let ((value (read-string "Value of pawproj: ")))
    (insert "pawproj = " value))
  )


(defun projwfc-plotboxes ()
  (interactive)
  (let ((value (read-string "Value of plotboxes: ")))
    (insert "plotboxes = " value))
  )


(defun projwfc-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun projwfc-tdosinboxes ()
  (interactive)
  (let ((value (read-string "Value of tdosinboxes: ")))
    (insert "tdosinboxes = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; dos- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun dos-DOS ()
  (interactive)
  (insert "&DOS")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; dos- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun dos-bz_sum ()
  (interactive)
  (let ((value (read-string "Value of bz_sum: ")))
    (insert "bz_sum = '" value "'"))
  (backward-char 1)
  )


(defun dos-degauss ()
  (interactive)
  (let ((value (read-string "Value of degauss: ")))
    (insert "degauss = " value))
  )


(defun dos-deltae ()
  (interactive)
  (let ((value (read-string "Value of DeltaE: ")))
    (insert "DeltaE = " value))
  )


(defun dos-emax ()
  (interactive)
  (let ((value (read-string "Value of Emax: ")))
    (insert "Emax = " value))
  )


(defun dos-emin ()
  (interactive)
  (let ((value (read-string "Value of Emin: ")))
    (insert "Emin = " value))
  )


(defun dos-fildos ()
  (interactive)
  (let ((value (read-string "Value of fildos: ")))
    (insert "fildos = '" value "'"))
  (backward-char 1)
  )


(defun dos-ngauss ()
  (interactive)
  (let ((value (read-string "Value of ngauss: ")))
    (insert "ngauss = " value))
  )


(defun dos-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun dos-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; bands- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun bands-BANDS ()
  (interactive)
  (insert "&BANDS")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; bands- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun bands-filband ()
  (interactive)
  (let ((value (read-string "Value of filband: ")))
    (insert "filband = '" value "'"))
  (backward-char 1)
  )


(defun bands-filp ()
  (interactive)
  (let ((value (read-string "Value of filp: ")))
    (insert "filp = '" value "'"))
  (backward-char 1)
  )


(defun bands-firstk ()
  (interactive)
  (let ((value (read-string "Value of firstk: ")))
    (insert "firstk = " value))
  )


(defun bands-lastk ()
  (interactive)
  (let ((value (read-string "Value of lastk: ")))
    (insert "lastk = " value))
  )


(defun bands-lp ()
  (interactive)
  (let ((value (read-string "Value of lp: ")))
    (insert "lp = " value))
  )


(defun bands-lsigma ()
  (interactive)
  (let ((value (read-string "Value of lsigma: ")))
    (insert "lsigma = " value))
  )


(defun bands-lsym ()
  (interactive)
  (let ((value (read-string "Value of lsym: ")))
    (insert "lsym = " value))
  )


(defun bands-no_overlap ()
  (interactive)
  (let ((value (read-string "Value of no_overlap: ")))
    (insert "no_overlap = " value))
  )


(defun bands-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun bands-plot_2d ()
  (interactive)
  (let ((value (read-string "Value of plot_2d: ")))
    (insert "plot_2d = " value))
  )


(defun bands-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun bands-spin_component ()
  (interactive)
  (let ((value (read-string "Value of spin_component: ")))
    (insert "spin_component = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pw2gw- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pw2gw-INPUTPP ()
  (interactive)
  (insert "&INPUTPP")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pw2gw- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pw2gw-deltae ()
  (interactive)
  (let ((value (read-string "Value of DeltaE: ")))
    (insert "DeltaE = " value))
  )


(defun pw2gw-emax ()
  (interactive)
  (let ((value (read-string "Value of Emax: ")))
    (insert "Emax = " value))
  )


(defun pw2gw-emin ()
  (interactive)
  (let ((value (read-string "Value of Emin: ")))
    (insert "Emin = " value))
  )


(defun pw2gw-outdir ()
  (interactive)
  (let ((value (read-string "Value of outdir: ")))
    (insert "outdir = " value))
  )


(defun pw2gw-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = " value))
  )


(defun pw2gw-qplda ()
  (interactive)
  (let ((value (read-string "Value of qplda: ")))
    (insert "qplda = " value))
  )


(defun pw2gw-vkb ()
  (interactive)
  (let ((value (read-string "Value of vkb: ")))
    (insert "vkb = " value))
  )


(defun pw2gw-vxcdiag ()
  (interactive)
  (let ((value (read-string "Value of vxcdiag: ")))
    (insert "vxcdiag = " value))
  )


(defun pw2gw-what ()
  (interactive)
  (let ((value (read-string "Value of what: ")))
    (insert "what = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pw2bgw- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pw2bgw-INPUT_PW2BGW ()
  (interactive)
  (insert "&INPUT_PW2BGW")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pw2bgw- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pw2bgw-outdir ()
  (interactive)
  (let ((value (read-string "Value of outdir: ")))
    (insert "outdir = " value))
  )


(defun pw2bgw-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = " value))
  )


(defun pw2bgw-real_or_complex ()
  (interactive)
  (let ((value (read-string "Value of real_or_complex: ")))
    (insert "real_or_complex = " value))
  )


(defun pw2bgw-rhog_file ()
  (interactive)
  (let ((value (read-string "Value of rhog_file: ")))
    (insert "rhog_file = " value))
  )


(defun pw2bgw-rhog_flag ()
  (interactive)
  (let ((value (read-string "Value of rhog_flag: ")))
    (insert "rhog_flag = " value))
  )


(defun pw2bgw-rhog_nvmax ()
  (interactive)
  (let ((value (read-string "Value of rhog_nvmax: ")))
    (insert "rhog_nvmax = " value))
  )


(defun pw2bgw-rhog_nvmin ()
  (interactive)
  (let ((value (read-string "Value of rhog_nvmin: ")))
    (insert "rhog_nvmin = " value))
  )


(defun pw2bgw-symm_type ()
  (interactive)
  (let ((value (read-string "Value of symm_type: ")))
    (insert "symm_type = " value))
  )


(defun pw2bgw-vkbg_file ()
  (interactive)
  (let ((value (read-string "Value of vkbg_file: ")))
    (insert "vkbg_file = " value))
  )


(defun pw2bgw-vkbg_flag ()
  (interactive)
  (let ((value (read-string "Value of vkbg_flag: ")))
    (insert "vkbg_flag = " value))
  )


(defun pw2bgw-vscg_file ()
  (interactive)
  (let ((value (read-string "Value of vscg_file: ")))
    (insert "vscg_file = " value))
  )


(defun pw2bgw-vscg_flag ()
  (interactive)
  (let ((value (read-string "Value of vscg_flag: ")))
    (insert "vscg_flag = " value))
  )


(defun pw2bgw-vxc0_file ()
  (interactive)
  (let ((value (read-string "Value of vxc0_file: ")))
    (insert "vxc0_file = " value))
  )


(defun pw2bgw-vxc0_flag ()
  (interactive)
  (let ((value (read-string "Value of vxc0_flag: ")))
    (insert "vxc0_flag = " value))
  )


(defun pw2bgw-vxc_diag_nmax ()
  (interactive)
  (let ((value (read-string "Value of vxc_diag_nmax: ")))
    (insert "vxc_diag_nmax = " value))
  )


(defun pw2bgw-vxc_diag_nmin ()
  (interactive)
  (let ((value (read-string "Value of vxc_diag_nmin: ")))
    (insert "vxc_diag_nmin = " value))
  )


(defun pw2bgw-vxc_file ()
  (interactive)
  (let ((value (read-string "Value of vxc_file: ")))
    (insert "vxc_file = " value))
  )


(defun pw2bgw-vxc_flag ()
  (interactive)
  (let ((value (read-string "Value of vxc_flag: ")))
    (insert "vxc_flag = " value))
  )


(defun pw2bgw-vxc_integral ()
  (interactive)
  (let ((value (read-string "Value of vxc_integral: ")))
    (insert "vxc_integral = " value))
  )


(defun pw2bgw-vxc_offdiag_nmax ()
  (interactive)
  (let ((value (read-string "Value of vxc_offdiag_nmax: ")))
    (insert "vxc_offdiag_nmax = " value))
  )


(defun pw2bgw-vxc_offdiag_nmin ()
  (interactive)
  (let ((value (read-string "Value of vxc_offdiag_nmin: ")))
    (insert "vxc_offdiag_nmin = " value))
  )


(defun pw2bgw-vxc_zero_rho_core ()
  (interactive)
  (let ((value (read-string "Value of vxc_zero_rho_core: ")))
    (insert "vxc_zero_rho_core = " value))
  )


(defun pw2bgw-vxcg_file ()
  (interactive)
  (let ((value (read-string "Value of vxcg_file: ")))
    (insert "vxcg_file = " value))
  )


(defun pw2bgw-vxcg_flag ()
  (interactive)
  (let ((value (read-string "Value of vxcg_flag: ")))
    (insert "vxcg_flag = " value))
  )


(defun pw2bgw-wfng_dk1 ()
  (interactive)
  (let ((value (read-string "Value of wfng_dk1: ")))
    (insert "wfng_dk1 = " value))
  )


(defun pw2bgw-wfng_dk2 ()
  (interactive)
  (let ((value (read-string "Value of wfng_dk2: ")))
    (insert "wfng_dk2 = " value))
  )


(defun pw2bgw-wfng_dk3 ()
  (interactive)
  (let ((value (read-string "Value of wfng_dk3: ")))
    (insert "wfng_dk3 = " value))
  )


(defun pw2bgw-wfng_file ()
  (interactive)
  (let ((value (read-string "Value of wfng_file: ")))
    (insert "wfng_file = " value))
  )


(defun pw2bgw-wfng_flag ()
  (interactive)
  (let ((value (read-string "Value of wfng_flag: ")))
    (insert "wfng_flag = " value))
  )


(defun pw2bgw-wfng_kgrid ()
  (interactive)
  (let ((value (read-string "Value of wfng_kgrid: ")))
    (insert "wfng_kgrid = " value))
  )


(defun pw2bgw-wfng_nk1 ()
  (interactive)
  (let ((value (read-string "Value of wfng_nk1: ")))
    (insert "wfng_nk1 = " value))
  )


(defun pw2bgw-wfng_nk2 ()
  (interactive)
  (let ((value (read-string "Value of wfng_nk2: ")))
    (insert "wfng_nk2 = " value))
  )


(defun pw2bgw-wfng_nk3 ()
  (interactive)
  (let ((value (read-string "Value of wfng_nk3: ")))
    (insert "wfng_nk3 = " value))
  )


(defun pw2bgw-wfng_nvmax ()
  (interactive)
  (let ((value (read-string "Value of wfng_nvmax: ")))
    (insert "wfng_nvmax = " value))
  )


(defun pw2bgw-wfng_nvmin ()
  (interactive)
  (let ((value (read-string "Value of wfng_nvmin: ")))
    (insert "wfng_nvmin = " value))
  )


(defun pw2bgw-wfng_occupation ()
  (interactive)
  (let ((value (read-string "Value of wfng_occupation: ")))
    (insert "wfng_occupation = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; molecularpdos- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun molecularpdos-INPUTMOPDOS ()
  (interactive)
  (insert "&INPUTMOPDOS")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; molecularpdos- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun molecularpdos-degauss ()
  (interactive)
  (let ((value (read-string "Value of degauss: ")))
    (insert "degauss = " value))
  )


(defun molecularpdos-deltae ()
  (interactive)
  (let ((value (read-string "Value of DeltaE: ")))
    (insert "DeltaE = " value))
  )


(defun molecularpdos-emax ()
  (interactive)
  (let ((value (read-string "Value of Emax: ")))
    (insert "Emax = " value))
  )


(defun molecularpdos-emin ()
  (interactive)
  (let ((value (read-string "Value of Emin: ")))
    (insert "Emin = " value))
  )


(defun molecularpdos-fileout ()
  (interactive)
  (let ((value (read-file-name "Value of fileout: ")))
    (insert "fileout = '" value "'"))
  (backward-char 1)
  )


(defun molecularpdos-i_atmwfc_beg_full ()
  (interactive)
  (let ((value (read-string "Value of i_atmwfc_beg_full: ")))
    (insert "i_atmwfc_beg_full = " value))
  )


(defun molecularpdos-i_atmwfc_beg_part ()
  (interactive)
  (let ((value (read-string "Value of i_atmwfc_beg_part: ")))
    (insert "i_atmwfc_beg_part = " value))
  )


(defun molecularpdos-i_atmwfc_end_full ()
  (interactive)
  (let ((value (read-string "Value of i_atmwfc_end_full: ")))
    (insert "i_atmwfc_end_full = " value))
  )


(defun molecularpdos-i_atmwfc_end_part ()
  (interactive)
  (let ((value (read-string "Value of i_atmwfc_end_part: ")))
    (insert "i_atmwfc_end_part = " value))
  )


(defun molecularpdos-i_bnd_beg_full ()
  (interactive)
  (let ((value (read-string "Value of i_bnd_beg_full: ")))
    (insert "i_bnd_beg_full = " value))
  )


(defun molecularpdos-i_bnd_beg_part ()
  (interactive)
  (let ((value (read-string "Value of i_bnd_beg_part: ")))
    (insert "i_bnd_beg_part = " value))
  )


(defun molecularpdos-i_bnd_end_full ()
  (interactive)
  (let ((value (read-string "Value of i_bnd_end_full: ")))
    (insert "i_bnd_end_full = " value))
  )


(defun molecularpdos-i_bnd_end_part ()
  (interactive)
  (let ((value (read-string "Value of i_bnd_end_part: ")))
    (insert "i_bnd_end_part = " value))
  )


(defun molecularpdos-kresolveddos ()
  (interactive)
  (let ((value (read-string "Value of kresolveddos: ")))
    (insert "kresolveddos = " value))
  )


(defun molecularpdos-ngauss ()
  (interactive)
  (let ((value (read-string "Value of ngauss: ")))
    (insert "ngauss = " value))
  )


(defun molecularpdos-xmlfile_full ()
  (interactive)
  (let ((value (read-string "Value of xmlfile_full: ")))
    (insert "xmlfile_full = " value))
  )


(defun molecularpdos-xmlfile_part ()
  (interactive)
  (let ((value (read-string "Value of xmlfile_part: ")))
    (insert "xmlfile_part = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pw_export- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pw_export-INPUTPP ()
  (interactive)
  (insert "&INPUTPP")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pw_export- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pw_export-ascii ()
  (interactive)
  (let ((value (read-string "Value of ascii: ")))
    (insert "ascii = " value))
  )


(defun pw_export-outdir ()
  (interactive)
  (let ((value (read-string "Value of outdir: ")))
    (insert "outdir = " value))
  )


(defun pw_export-pp_file ()
  (interactive)
  (let ((value (read-string "Value of pp_file: ")))
    (insert "pp_file = " value))
  )


(defun pw_export-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = " value))
  )


(defun pw_export-pseudo_dir ()
  (interactive)
  (let ((value (read-string "Value of pseudo_dir: ")))
    (insert "pseudo_dir = " value))
  )


(defun pw_export-psfile ()
  (interactive)
  (let ((value (read-string "Value of psfile: ")))
    (insert "psfile = " value))
  )


(defun pw_export-single_file ()
  (interactive)
  (let ((value (read-string "Value of single_file: ")))
    (insert "single_file = " value))
  )


(defun pw_export-uspp_spsi ()
  (interactive)
  (let ((value (read-string "Value of uspp_spsi: ")))
    (insert "uspp_spsi = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; bgw2pw- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun bgw2pw-INPUT_BGW2PW ()
  (interactive)
  (insert "&INPUT_BGW2PW")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; bgw2pw- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun bgw2pw-outdir ()
  (interactive)
  (let ((value (read-string "Value of outdir: ")))
    (insert "outdir = " value))
  )


(defun bgw2pw-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = " value))
  )


(defun bgw2pw-real_or_complex ()
  (interactive)
  (let ((value (read-string "Value of real_or_complex: ")))
    (insert "real_or_complex = " value))
  )


(defun bgw2pw-rhog_file ()
  (interactive)
  (let ((value (read-string "Value of rhog_file: ")))
    (insert "rhog_file = " value))
  )


(defun bgw2pw-rhog_flag ()
  (interactive)
  (let ((value (read-string "Value of rhog_flag: ")))
    (insert "rhog_flag = " value))
  )


(defun bgw2pw-wfng_file ()
  (interactive)
  (let ((value (read-string "Value of wfng_file: ")))
    (insert "wfng_file = " value))
  )


(defun bgw2pw-wfng_flag ()
  (interactive)
  (let ((value (read-string "Value of wfng_flag: ")))
    (insert "wfng_flag = " value))
  )


(defun bgw2pw-wfng_nband ()
  (interactive)
  (let ((value (read-string "Value of wfng_nband: ")))
    (insert "wfng_nband = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pp- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pp-INPUTPP ()
  (interactive)
  (insert "&INPUTPP")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun pp-PLOT ()
  (interactive)
  (insert "&PLOT")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pp- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pp-degauss_ldos ()
  (interactive)
  (let ((value (read-string "Value of degauss_ldos: ")))
    (insert "degauss_ldos = " value))
  )


(defun pp-delta_e ()
  (interactive)
  (let ((value (read-string "Value of delta_e: ")))
    (insert "delta_e = " value))
  )


(defun pp-e1 ()
  (interactive)
  (let ((value (read-string "Value of e1: ")))
    (insert "e1 = " value))
  )


(defun pp-e2 ()
  (interactive)
  (let ((value (read-string "Value of e2: ")))
    (insert "e2 = " value))
  )


(defun pp-e3 ()
  (interactive)
  (let ((value (read-string "Value of e3: ")))
    (insert "e3 = " value))
  )


(defun pp-emax ()
  (interactive)
  (let ((value (read-string "Value of emax: ")))
    (insert "emax = " value))
  )


(defun pp-emin ()
  (interactive)
  (let ((value (read-string "Value of emin: ")))
    (insert "emin = " value))
  )


(defun pp-fileout ()
  (interactive)
  (let ((value (read-file-name "Value of fileout: ")))
    (insert "fileout = '" value "'"))
  (backward-char 1)
  )


(defun pp-filepp ()
  (interactive)
  (let ((value (read-file-name "Value of filepp: ")))
    (insert "filepp = '" value "'"))
  (backward-char 1)
  )


(defun pp-filplot ()
  (interactive)
  (let ((value (read-string "Value of filplot: ")))
    (insert "filplot = '" value "'"))
  (backward-char 1)
  )


(defun pp-iflag ()
  (interactive)
  (let ((value (read-string "Value of iflag: ")))
    (insert "iflag = " value))
  )


(defun pp-interpolation ()
  (interactive)
  (let ((value (read-string "Value of interpolation: ")))
    (insert "interpolation = '" value "'"))
  (backward-char 1)
  )


(defun pp-kband ()
  (interactive)
  (let ((value (read-string "Value of kband: ")))
    (insert "kband = " value))
  )


(defun pp-kpoint ()
  (interactive)
  (let ((value (read-string "Value of kpoint: ")))
    (insert "kpoint = " value))
  )


(defun pp-lsign ()
  (interactive)
  (let ((value (read-string "Value of lsign: ")))
    (insert "lsign = " value))
  )


(defun pp-nfile ()
  (interactive)
  (let ((value (read-string "Value of nfile: ")))
    (insert "nfile = " value))
  )


(defun pp-nx ()
  (interactive)
  (let ((value (read-string "Value of nx: ")))
    (insert "nx = " value))
  )


(defun pp-ny ()
  (interactive)
  (let ((value (read-string "Value of ny: ")))
    (insert "ny = " value))
  )


(defun pp-nz ()
  (interactive)
  (let ((value (read-string "Value of nz: ")))
    (insert "nz = " value))
  )


(defun pp-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun pp-output_format ()
  (interactive)
  (let ((value (read-string "Value of output_format: ")))
    (insert "output_format = " value))
  )


(defun pp-plot_num ()
  (interactive)
  (let ((value (read-string "Value of plot_num: ")))
    (insert "plot_num = " value))
  )


(defun pp-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun pp-radius ()
  (interactive)
  (let ((value (read-string "Value of radius: ")))
    (insert "radius = " value))
  )


(defun pp-sample_bias ()
  (interactive)
  (let ((value (read-string "Value of sample_bias: ")))
    (insert "sample_bias = " value))
  )


(defun pp-spin_component ()
  (interactive)
  (let ((value (read-string "Value of spin_component: ")))
    (insert "spin_component = " value))
  )


(defun pp-weight ()
  (interactive)
  (let ((value (read-string "Value of weight: ")))
    (insert "weight = " value))
  )


(defun pp-x0 ()
  (interactive)
  (let ((value (read-string "Value of x0: ")))
    (insert "x0 = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; hp- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun hp-INPUTHP ()
  (interactive)
  (insert "&INPUTHP")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; hp- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun hp-alpha_mix ()
  (interactive)
  (let ((value (read-string "Value of alpha_mix: ")))
    (insert "alpha_mix = " value))
  )


(defun hp-compute_hp ()
  (interactive)
  (let ((value (read-string "Value of compute_hp: ")))
    (insert "compute_hp = " value))
  )


(defun hp-conv_thr_chi ()
  (interactive)
  (let ((value (read-string "Value of conv_thr_chi: ")))
    (insert "conv_thr_chi = " value))
  )


(defun hp-determine_num_pert_only ()
  (interactive)
  (let ((value (read-string "Value of determine_num_pert_only: ")))
    (insert "determine_num_pert_only = " value))
  )


(defun hp-docc_thr ()
  (interactive)
  (let ((value (read-string "Value of docc_thr: ")))
    (insert "docc_thr = " value))
  )


(defun hp-equiv_type ()
  (interactive)
  (let ((value (read-string "Value of equiv_type: ")))
    (insert "equiv_type = " value))
  )


(defun hp-ethr_nscf ()
  (interactive)
  (let ((value (read-string "Value of ethr_nscf: ")))
    (insert "ethr_nscf = " value))
  )


(defun hp-find_atpert ()
  (interactive)
  (let ((value (read-string "Value of find_atpert: ")))
    (insert "find_atpert = " value))
  )


(defun hp-iverbosity ()
  (interactive)
  (let ((value (read-string "Value of iverbosity: ")))
    (insert "iverbosity = " value))
  )


(defun hp-last_q ()
  (interactive)
  (let ((value (read-string "Value of last_q: ")))
    (insert "last_q = " value))
  )


(defun hp-max_seconds ()
  (interactive)
  (let ((value (read-string "Value of max_seconds: ")))
    (insert "max_seconds = " value))
  )


(defun hp-niter_max ()
  (interactive)
  (let ((value (read-string "Value of niter_max: ")))
    (insert "niter_max = " value))
  )


(defun hp-nmix ()
  (interactive)
  (let ((value (read-string "Value of nmix: ")))
    (insert "nmix = " value))
  )


(defun hp-nq1 ()
  (interactive)
  (let ((value (read-string "Value of nq1: ")))
    (insert "nq1 = " value))
  )


(defun hp-nq2 ()
  (interactive)
  (let ((value (read-string "Value of nq2: ")))
    (insert "nq2 = " value))
  )


(defun hp-nq3 ()
  (interactive)
  (let ((value (read-string "Value of nq3: ")))
    (insert "nq3 = " value))
  )


(defun hp-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun hp-perturb_only_atom ()
  (interactive)
  (let ((value (read-string "Value of perturb_only_atom: ")))
    (insert "perturb_only_atom = " value))
  )


(defun hp-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun hp-skip_equivalence_q ()
  (interactive)
  (let ((value (read-string "Value of skip_equivalence_q: ")))
    (insert "skip_equivalence_q = " value))
  )


(defun hp-skip_type ()
  (interactive)
  (let ((value (read-string "Value of skip_type: ")))
    (insert "skip_type = " value))
  )


(defun hp-start_q ()
  (interactive)
  (let ((value (read-string "Value of start_q: ")))
    (insert "start_q = " value))
  )


(defun hp-sum_pertq ()
  (interactive)
  (let ((value (read-string "Value of sum_pertq: ")))
    (insert "sum_pertq = " value))
  )


(defun hp-thresh_init ()
  (interactive)
  (let ((value (read-string "Value of thresh_init: ")))
    (insert "thresh_init = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; ld1- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ld1-INPUT ()
  (interactive)
  (insert "&INPUT")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun ld1-INPUTP ()
  (interactive)
  (insert "&INPUTP")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun ld1-TEST ()
  (interactive)
  (insert "&TEST")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; ld1- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ld1-atom ()
  (interactive)
  (let ((value (read-string "Value of atom: ")))
    (insert "atom = '" value "'"))
  (backward-char 1)
  )


(defun ld1-author ()
  (interactive)
  (let ((value (read-string "Value of author: ")))
    (insert "author = '" value "'"))
  (backward-char 1)
  )


(defun ld1-beta ()
  (interactive)
  (let ((value (read-string "Value of beta: ")))
    (insert "beta = " value))
  )


(defun ld1-cau_fact ()
  (interactive)
  (let ((value (read-string "Value of cau_fact: ")))
    (insert "cau_fact = " value))
  )


(defun ld1-config ()
  (interactive)
  (let ((value (read-string "Value of config: ")))
    (insert "config = '" value "'"))
  (backward-char 1)
  )


(defun ld1-configts ()
  (interactive)
  (let ((value (read-string "Value of configts: ")))
    (insert "configts = '" value "'"))
  (backward-char 1)
  )


(defun ld1-decut ()
  (interactive)
  (let ((value (read-string "Value of decut: ")))
    (insert "decut = " value))
  )


(defun ld1-deld ()
  (interactive)
  (let ((value (read-string "Value of deld: ")))
    (insert "deld = " value))
  )


(defun ld1-dft ()
  (interactive)
  (let ((value (read-string "Value of dft: ")))
    (insert "dft = '" value "'"))
  (backward-char 1)
  )


(defun ld1-dx ()
  (interactive)
  (let ((value (read-string "Value of dx: ")))
    (insert "dx = " value))
  )


(defun ld1-ecutmax ()
  (interactive)
  (let ((value (read-string "Value of ecutmax: ")))
    (insert "ecutmax = " value))
  )


(defun ld1-ecutmin ()
  (interactive)
  (let ((value (read-string "Value of ecutmin: ")))
    (insert "ecutmin = " value))
  )


(defun ld1-emaxld ()
  (interactive)
  (let ((value (read-string "Value of emaxld: ")))
    (insert "emaxld = " value))
  )


(defun ld1-eminld ()
  (interactive)
  (let ((value (read-string "Value of eminld: ")))
    (insert "eminld = " value))
  )


(defun ld1-file_beta ()
  (interactive)
  (let ((value (read-file-name "Value of file_beta: ")))
    (insert "file_beta = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_charge ()
  (interactive)
  (let ((value (read-file-name "Value of file_charge: ")))
    (insert "file_charge = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_chi ()
  (interactive)
  (let ((value (read-file-name "Value of file_chi: ")))
    (insert "file_chi = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_core ()
  (interactive)
  (let ((value (read-file-name "Value of file_core: ")))
    (insert "file_core = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_pseudo ()
  (interactive)
  (let ((value (read-file-name "Value of file_pseudo: ")))
    (insert "file_pseudo = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_pseudopw ()
  (interactive)
  (let ((value (read-file-name "Value of file_pseudopw: ")))
    (insert "file_pseudopw = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_qvan ()
  (interactive)
  (let ((value (read-file-name "Value of file_qvan: ")))
    (insert "file_qvan = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_recon ()
  (interactive)
  (let ((value (read-file-name "Value of file_recon: ")))
    (insert "file_recon = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_screen ()
  (interactive)
  (let ((value (read-file-name "Value of file_screen: ")))
    (insert "file_screen = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_wfcaegen ()
  (interactive)
  (let ((value (read-file-name "Value of file_wfcaegen: ")))
    (insert "file_wfcaegen = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_wfcncgen ()
  (interactive)
  (let ((value (read-file-name "Value of file_wfcncgen: ")))
    (insert "file_wfcncgen = '" value "'"))
  (backward-char 1)
  )


(defun ld1-file_wfcusgen ()
  (interactive)
  (let ((value (read-file-name "Value of file_wfcusgen: ")))
    (insert "file_wfcusgen = '" value "'"))
  (backward-char 1)
  )


(defun ld1-frozen_core ()
  (interactive)
  (let ((value (read-string "Value of frozen_core: ")))
    (insert "frozen_core = " value))
  )


(defun ld1-isic ()
  (interactive)
  (let ((value (read-string "Value of isic: ")))
    (insert "isic = " value))
  )


(defun ld1-iswitch ()
  (interactive)
  (let ((value (read-string "Value of iswitch: ")))
    (insert "iswitch = " value))
  )


(defun ld1-latt ()
  (interactive)
  (let ((value (read-string "Value of latt: ")))
    (insert "latt = " value))
  )


(defun ld1-lgipaw_reconstruction ()
  (interactive)
  (let ((value (read-string "Value of lgipaw_reconstruction: ")))
    (insert "lgipaw_reconstruction = " value))
  )


(defun ld1-lloc ()
  (interactive)
  (let ((value (read-string "Value of lloc: ")))
    (insert "lloc = " value))
  )


(defun ld1-lpaw ()
  (interactive)
  (let ((value (read-string "Value of lpaw: ")))
    (insert "lpaw = " value))
  )


(defun ld1-lsave_wfc ()
  (interactive)
  (let ((value (read-string "Value of lsave_wfc: ")))
    (insert "lsave_wfc = " value))
  )


(defun ld1-lsd ()
  (interactive)
  (let ((value (read-string "Value of lsd: ")))
    (insert "lsd = " value))
  )


(defun ld1-lsdts ()
  (interactive)
  (let ((value (read-string "Value of lsdts: ")))
    (insert "lsdts = " value))
  )


(defun ld1-lsmall ()
  (interactive)
  (let ((value (read-string "Value of lsmall: ")))
    (insert "lsmall = " value))
  )


(defun ld1-max_out_wfc ()
  (interactive)
  (let ((value (read-string "Value of max_out_wfc: ")))
    (insert "max_out_wfc = " value))
  )


(defun ld1-nconf ()
  (interactive)
  (let ((value (read-string "Value of nconf: ")))
    (insert "nconf = " value))
  )


(defun ld1-new_core_ps ()
  (interactive)
  (let ((value (read-string "Value of new_core_ps: ")))
    (insert "new_core_ps = " value))
  )


(defun ld1-nlcc ()
  (interactive)
  (let ((value (read-string "Value of nlcc: ")))
    (insert "nlcc = " value))
  )


(defun ld1-nld ()
  (interactive)
  (let ((value (read-string "Value of nld: ")))
    (insert "nld = " value))
  )


(defun ld1-noscf ()
  (interactive)
  (let ((value (read-string "Value of noscf: ")))
    (insert "noscf = " value))
  )


(defun ld1-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun ld1-pseudotype ()
  (interactive)
  (let ((value (read-string "Value of pseudotype: ")))
    (insert "pseudotype = " value))
  )


(defun ld1-rcloc ()
  (interactive)
  (let ((value (read-string "Value of rcloc: ")))
    (insert "rcloc = " value))
  )


(defun ld1-rcore ()
  (interactive)
  (let ((value (read-string "Value of rcore: ")))
    (insert "rcore = " value))
  )


(defun ld1-rcutv ()
  (interactive)
  (let ((value (read-string "Value of rcutv: ")))
    (insert "rcutv = " value))
  )


(defun ld1-rel ()
  (interactive)
  (let ((value (read-string "Value of rel: ")))
    (insert "rel = " value))
  )


(defun ld1-rel_dist ()
  (interactive)
  (let ((value (read-string "Value of rel_dist: ")))
    (insert "rel_dist = '" value "'"))
  (backward-char 1)
  )


(defun ld1-relpert ()
  (interactive)
  (let ((value (read-string "Value of relpert: ")))
    (insert "relpert = " value))
  )


(defun ld1-rho0 ()
  (interactive)
  (let ((value (read-string "Value of rho0: ")))
    (insert "rho0 = " value))
  )


(defun ld1-rlderiv ()
  (interactive)
  (let ((value (read-string "Value of rlderiv: ")))
    (insert "rlderiv = " value))
  )


(defun ld1-rm ()
  (interactive)
  (let ((value (read-string "Value of rm: ")))
    (insert "rm = " value))
  )


(defun ld1-rmatch_augfun ()
  (interactive)
  (let ((value (read-string "Value of rmatch_augfun: ")))
    (insert "rmatch_augfun = " value))
  )


(defun ld1-rmatch_augfun_nc ()
  (interactive)
  (let ((value (read-string "Value of rmatch_augfun_nc: ")))
    (insert "rmatch_augfun_nc = " value))
  )


(defun ld1-rmax ()
  (interactive)
  (let ((value (read-string "Value of rmax: ")))
    (insert "rmax = " value))
  )


(defun ld1-rpwe ()
  (interactive)
  (let ((value (read-string "Value of rpwe: ")))
    (insert "rpwe = " value))
  )


(defun ld1-rytoev_fact ()
  (interactive)
  (let ((value (read-string "Value of rytoev_fact: ")))
    (insert "rytoev_fact = " value))
  )


(defun ld1-title ()
  (interactive)
  (let ((value (read-string "Value of title: ")))
    (insert "title = '" value "'"))
  (backward-char 1)
  )


(defun ld1-tm ()
  (interactive)
  (let ((value (read-string "Value of tm: ")))
    (insert "tm = " value))
  )


(defun ld1-tr2 ()
  (interactive)
  (let ((value (read-string "Value of tr2: ")))
    (insert "tr2 = " value))
  )


(defun ld1-use_paw_as_gipaw ()
  (interactive)
  (let ((value (read-string "Value of use_paw_as_gipaw: ")))
    (insert "use_paw_as_gipaw = " value))
  )


(defun ld1-vdw ()
  (interactive)
  (let ((value (read-string "Value of vdw: ")))
    (insert "vdw = " value))
  )


(defun ld1-verbosity ()
  (interactive)
  (let ((value (read-string "Value of verbosity: ")))
    (insert "verbosity = '" value "'"))
  (backward-char 1)
  )


(defun ld1-which_augfun ()
  (interactive)
  (let ((value (read-string "Value of which_augfun: ")))
    (insert "which_augfun = '" value "'"))
  (backward-char 1)
  )


(defun ld1-write_coulomb ()
  (interactive)
  (let ((value (read-string "Value of write_coulomb: ")))
    (insert "write_coulomb = " value))
  )


(defun ld1-xmin ()
  (interactive)
  (let ((value (read-string "Value of xmin: ")))
    (insert "xmin = " value))
  )


(defun ld1-zed ()
  (interactive)
  (let ((value (read-string "Value of zed: ")))
    (insert "zed = " value))
  )


(defun ld1-zval ()
  (interactive)
  (let ((value (read-string "Value of zval: ")))
    (insert "zval = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; neb- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun neb-PATH ()
  (interactive)
  (insert "&PATH")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; neb- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun neb-ci_scheme ()
  (interactive)
  (let ((value (read-string "Value of CI_scheme: ")))
    (insert "CI_scheme = '" value "'"))
  (backward-char 1)
  )


(defun neb-ds ()
  (interactive)
  (let ((value (read-string "Value of ds: ")))
    (insert "ds = " value))
  )


(defun neb-fcp_mu ()
  (interactive)
  (let ((value (read-string "Value of fcp_mu: ")))
    (insert "fcp_mu = " value))
  )


(defun neb-fcp_tot_charge_first ()
  (interactive)
  (let ((value (read-string "Value of fcp_tot_charge_first: ")))
    (insert "fcp_tot_charge_first = " value))
  )


(defun neb-fcp_tot_charge_last ()
  (interactive)
  (let ((value (read-string "Value of fcp_tot_charge_last: ")))
    (insert "fcp_tot_charge_last = " value))
  )


(defun neb-first_last_opt ()
  (interactive)
  (let ((value (read-string "Value of first_last_opt: ")))
    (insert "first_last_opt = " value))
  )


(defun neb-k_max ()
  (interactive)
  (let ((value (read-string "Value of k_max: ")))
    (insert "k_max = " value))
  )


(defun neb-k_min ()
  (interactive)
  (let ((value (read-string "Value of k_min: ")))
    (insert "k_min = " value))
  )


(defun neb-lfcpopt ()
  (interactive)
  (let ((value (read-string "Value of lfcpopt: ")))
    (insert "lfcpopt = " value))
  )


(defun neb-minimum_image ()
  (interactive)
  (let ((value (read-string "Value of minimum_image: ")))
    (insert "minimum_image = " value))
  )


(defun neb-nstep_path ()
  (interactive)
  (let ((value (read-string "Value of nstep_path: ")))
    (insert "nstep_path = " value))
  )


(defun neb-num_of_images ()
  (interactive)
  (let ((value (read-string "Value of num_of_images: ")))
    (insert "num_of_images = " value))
  )


(defun neb-opt_scheme ()
  (interactive)
  (let ((value (read-string "Value of opt_scheme: ")))
    (insert "opt_scheme = '" value "'"))
  (backward-char 1)
  )


(defun neb-path_thr ()
  (interactive)
  (let ((value (read-string "Value of path_thr: ")))
    (insert "path_thr = " value))
  )


(defun neb-restart_mode ()
  (interactive)
  (let ((value (read-string "Value of restart_mode: ")))
    (insert "restart_mode = '" value "'"))
  (backward-char 1)
  )


(defun neb-string_method ()
  (interactive)
  (let ((value (read-string "Value of string_method: ")))
    (insert "string_method = '" value "'"))
  (backward-char 1)
  )


(defun neb-temp_req ()
  (interactive)
  (let ((value (read-string "Value of temp_req: ")))
    (insert "temp_req = " value))
  )


(defun neb-use_freezing ()
  (interactive)
  (let ((value (read-string "Value of use_freezing: ")))
    (insert "use_freezing = " value))
  )


(defun neb-use_masses ()
  (interactive)
  (let ((value (read-string "Value of use_masses: ")))
    (insert "use_masses = " value))
  )


	   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; neb- supercards functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun neb-BEGIN ()
  (interactive)
  (insert "BEGIN")
  (newline 2)
  (insert "END")
  (forward-line -1)
  )


(defun neb-BEGIN_ENGINE_INPUT ()
  (interactive)
  (insert "BEGIN_ENGINE_INPUT")
  (newline 2)
  (insert "END_ENGINE_INPUT")
  (forward-line -1)
  )


(defun neb-BEGIN_PATH_INPUT ()
  (interactive)
  (insert "BEGIN_PATH_INPUT")
  (newline 2)
  (insert "END_PATH_INPUT")
  (forward-line -1)
  )


(defun neb-BEGIN_POSITIONS ()
  (interactive)
  (insert "BEGIN_POSITIONS")
  (newline 2)
  (insert "END_POSITIONS")
  (forward-line -1)
  )


(defun neb-FIRST_IMAGE ()
  (interactive)
  (insert "FIRST_IMAGE")
  (newline 1)
  )


(defun neb-INTERMEDIATE_IMAGE ()
  (interactive)
  (insert "INTERMEDIATE_IMAGE")
  (newline 1)
  )


(defun neb-LAST_IMAGE ()
  (interactive)
  (insert "LAST_IMAGE")
  (newline 1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; neb- cards functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun neb-ATOMIC_POSITIONS ()
  (interactive)
  (let ((flag (read-string "Flags: { alat | bohr | angstrom | crystal | crystal_sg } ")))
    (insert "ATOMIC_POSITIONS " flag))
  (newline 1)
  )


(defun neb-CLIMBING_IMAGES ()
 (interactive)
 (insert "CLIMBING_IMAGES")
 (newline 1)
 )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; ph- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ph-INPUTPH ()
  (interactive)
  (insert "&INPUTPH")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; ph- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun ph-alpha_mix ()
  (interactive)
  (let ((value (read-string "Value of alpha_mix: ")))
    (insert "alpha_mix = " value))
  )


(defun ph-amass ()
  (interactive)
  (let ((value (read-string "Value of amass: ")))
    (insert "amass = " value))
  )


(defun ph-asr ()
  (interactive)
  (let ((value (read-string "Value of asr: ")))
    (insert "asr = " value))
  )


(defun ph-dek ()
  (interactive)
  (let ((value (read-string "Value of dek: ")))
    (insert "dek = " value))
  )


(defun ph-drho_star ()
  (interactive)
  (let ((value (read-string "Value of drho_star: ")))
    (insert "drho_star = " value))
  )


(defun ph-dvscf_star ()
  (interactive)
  (let ((value (read-string "Value of dvscf_star: ")))
    (insert "dvscf_star = " value))
  )


(defun ph-electron_phonon ()
  (interactive)
  (let ((value (read-string "Value of electron_phonon: ")))
    (insert "electron_phonon = '" value "'"))
  (backward-char 1)
  )


(defun ph-elop ()
  (interactive)
  (let ((value (read-string "Value of elop: ")))
    (insert "elop = " value))
  )


(defun ph-epsil ()
  (interactive)
  (let ((value (read-string "Value of epsil: ")))
    (insert "epsil = " value))
  )


(defun ph-eth_ns ()
  (interactive)
  (let ((value (read-string "Value of eth_ns: ")))
    (insert "eth_ns = " value))
  )


(defun ph-eth_rps ()
  (interactive)
  (let ((value (read-string "Value of eth_rps: ")))
    (insert "eth_rps = " value))
  )


(defun ph-fildrho ()
  (interactive)
  (let ((value (read-string "Value of fildrho: ")))
    (insert "fildrho = '" value "'"))
  (backward-char 1)
  )


(defun ph-fildvscf ()
  (interactive)
  (let ((value (read-string "Value of fildvscf: ")))
    (insert "fildvscf = '" value "'"))
  (backward-char 1)
  )


(defun ph-fildyn ()
  (interactive)
  (let ((value (read-string "Value of fildyn: ")))
    (insert "fildyn = '" value "'"))
  (backward-char 1)
  )


(defun ph-fpol ()
  (interactive)
  (let ((value (read-string "Value of fpol: ")))
    (insert "fpol = " value))
  )


(defun ph-k1 ()
  (interactive)
  (let ((value (read-string "Value of k1: ")))
    (insert "k1 = " value))
  )


(defun ph-k2 ()
  (interactive)
  (let ((value (read-string "Value of k2: ")))
    (insert "k2 = " value))
  )


(defun ph-k3 ()
  (interactive)
  (let ((value (read-string "Value of k3: ")))
    (insert "k3 = " value))
  )


(defun ph-last_irr ()
  (interactive)
  (let ((value (read-string "Value of last_irr: ")))
    (insert "last_irr = " value))
  )


(defun ph-last_q ()
  (interactive)
  (let ((value (read-string "Value of last_q: ")))
    (insert "last_q = " value))
  )


(defun ph-ldiag ()
  (interactive)
  (let ((value (read-string "Value of ldiag: ")))
    (insert "ldiag = " value))
  )


(defun ph-ldisp ()
  (interactive)
  (let ((value (read-string "Value of ldisp: ")))
    (insert "ldisp = " value))
  )


(defun ph-lnoloc ()
  (interactive)
  (let ((value (read-string "Value of lnoloc: ")))
    (insert "lnoloc = " value))
  )


(defun ph-low_directory_check ()
  (interactive)
  (let ((value (read-string "Value of low_directory_check: ")))
    (insert "low_directory_check = " value))
  )


(defun ph-lqdir ()
  (interactive)
  (let ((value (read-string "Value of lqdir: ")))
    (insert "lqdir = " value))
  )


(defun ph-lraman ()
  (interactive)
  (let ((value (read-string "Value of lraman: ")))
    (insert "lraman = " value))
  )


(defun ph-lrpa ()
  (interactive)
  (let ((value (read-string "Value of lrpa: ")))
    (insert "lrpa = " value))
  )


(defun ph-lshift_q ()
  (interactive)
  (let ((value (read-string "Value of lshift_q: ")))
    (insert "lshift_q = " value))
  )


(defun ph-max_seconds ()
  (interactive)
  (let ((value (read-string "Value of max_seconds: ")))
    (insert "max_seconds = " value))
  )


(defun ph-modenum ()
  (interactive)
  (let ((value (read-string "Value of modenum: ")))
    (insert "modenum = " value))
  )


(defun ph-nat_todo ()
  (interactive)
  (let ((value (read-string "Value of nat_todo: ")))
    (insert "nat_todo = " value))
  )


(defun ph-niter_ph ()
  (interactive)
  (let ((value (read-string "Value of niter_ph: ")))
    (insert "niter_ph = " value))
  )


(defun ph-nk1 ()
  (interactive)
  (let ((value (read-string "Value of nk1: ")))
    (insert "nk1 = " value))
  )


(defun ph-nk2 ()
  (interactive)
  (let ((value (read-string "Value of nk2: ")))
    (insert "nk2 = " value))
  )


(defun ph-nk3 ()
  (interactive)
  (let ((value (read-string "Value of nk3: ")))
    (insert "nk3 = " value))
  )


(defun ph-nmix_ph ()
  (interactive)
  (let ((value (read-string "Value of nmix_ph: ")))
    (insert "nmix_ph = " value))
  )


(defun ph-nogg ()
  (interactive)
  (let ((value (read-string "Value of nogg: ")))
    (insert "nogg = " value))
  )


(defun ph-nq1 ()
  (interactive)
  (let ((value (read-string "Value of nq1: ")))
    (insert "nq1 = " value))
  )


(defun ph-nq2 ()
  (interactive)
  (let ((value (read-string "Value of nq2: ")))
    (insert "nq2 = " value))
  )


(defun ph-nq3 ()
  (interactive)
  (let ((value (read-string "Value of nq3: ")))
    (insert "nq3 = " value))
  )


(defun ph-only_init ()
  (interactive)
  (let ((value (read-string "Value of only_init: ")))
    (insert "only_init = " value))
  )


(defun ph-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun ph-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun ph-q2d ()
  (interactive)
  (let ((value (read-string "Value of q2d: ")))
    (insert "q2d = " value))
  )


(defun ph-q_in_band_form ()
  (interactive)
  (let ((value (read-string "Value of q_in_band_form: ")))
    (insert "q_in_band_form = " value))
  )


(defun ph-qplot ()
  (interactive)
  (let ((value (read-string "Value of qplot: ")))
    (insert "qplot = " value))
  )


(defun ph-read_dns_bare ()
  (interactive)
  (let ((value (read-string "Value of read_dns_bare: ")))
    (insert "read_dns_bare = " value))
  )


(defun ph-recover ()
  (interactive)
  (let ((value (read-string "Value of recover: ")))
    (insert "recover = " value))
  )


(defun ph-reduce_io ()
  (interactive)
  (let ((value (read-string "Value of reduce_io: ")))
    (insert "reduce_io = " value))
  )


(defun ph-search_sym ()
  (interactive)
  (let ((value (read-string "Value of search_sym: ")))
    (insert "search_sym = " value))
  )


(defun ph-start_irr ()
  (interactive)
  (let ((value (read-string "Value of start_irr: ")))
    (insert "start_irr = " value))
  )


(defun ph-start_q ()
  (interactive)
  (let ((value (read-string "Value of start_q: ")))
    (insert "start_q = " value))
  )


(defun ph-tr2_ph ()
  (interactive)
  (let ((value (read-string "Value of tr2_ph: ")))
    (insert "tr2_ph = " value))
  )


(defun ph-trans ()
  (interactive)
  (let ((value (read-string "Value of trans: ")))
    (insert "trans = " value))
  )


(defun ph-verbosity ()
  (interactive)
  (let ((value (read-string "Value of verbosity: ")))
    (insert "verbosity = '" value "'"))
  (backward-char 1)
  )


(defun ph-zeu ()
  (interactive)
  (let ((value (read-string "Value of zeu: ")))
    (insert "zeu = " value))
  )


(defun ph-zue ()
  (interactive)
  (let ((value (read-string "Value of zue: ")))
    (insert "zue = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; dynmat- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun dynmat-INPUT ()
  (interactive)
  (insert "&INPUT")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; dynmat- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun dynmat-amass ()
  (interactive)
  (let ((value (read-string "Value of amass: ")))
    (insert "amass = " value))
  )


(defun dynmat-asr ()
  (interactive)
  (let ((value (read-string "Value of asr: ")))
    (insert "asr = '" value "'"))
  (backward-char 1)
  )


(defun dynmat-axis ()
  (interactive)
  (let ((value (read-string "Value of axis: ")))
    (insert "axis = " value))
  )


(defun dynmat-fildyn ()
  (interactive)
  (let ((value (read-string "Value of fildyn: ")))
    (insert "fildyn = '" value "'"))
  (backward-char 1)
  )


(defun dynmat-fileig ()
  (interactive)
  (let ((value (read-file-name "Value of fileig: ")))
    (insert "fileig = '" value "'"))
  (backward-char 1)
  )


(defun dynmat-filmol ()
  (interactive)
  (let ((value (read-string "Value of filmol: ")))
    (insert "filmol = '" value "'"))
  (backward-char 1)
  )


(defun dynmat-filout ()
  (interactive)
  (let ((value (read-string "Value of filout: ")))
    (insert "filout = '" value "'"))
  (backward-char 1)
  )


(defun dynmat-filxsf ()
  (interactive)
  (let ((value (read-string "Value of filxsf: ")))
    (insert "filxsf = '" value "'"))
  (backward-char 1)
  )


(defun dynmat-loto_2d ()
  (interactive)
  (let ((value (read-string "Value of loto_2d: ")))
    (insert "loto_2d = " value))
  )


(defun dynmat-lperm ()
  (interactive)
  (let ((value (read-string "Value of lperm: ")))
    (insert "lperm = " value))
  )


(defun dynmat-lplasma ()
  (interactive)
  (let ((value (read-string "Value of lplasma: ")))
    (insert "lplasma = " value))
  )


(defun dynmat-q ()
  (interactive)
  (let ((value (read-string "Value of q: ")))
    (insert "q = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pwcond- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pwcond-INPUTCOND ()
  (interactive)
  (insert "&INPUTCOND")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; pwcond- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun pwcond-band_file ()
  (interactive)
  (let ((value (read-file-name "Value of band_file: ")))
    (insert "band_file = '" value "'"))
  (backward-char 1)
  )


(defun pwcond-bdl ()
  (interactive)
  (let ((value (read-string "Value of bdl: ")))
    (insert "bdl = " value))
  )


(defun pwcond-bdr ()
  (interactive)
  (let ((value (read-string "Value of bdr: ")))
    (insert "bdr = " value))
  )


(defun pwcond-bds ()
  (interactive)
  (let ((value (read-string "Value of bds: ")))
    (insert "bds = " value))
  )


(defun pwcond-denergy ()
  (interactive)
  (let ((value (read-string "Value of denergy: ")))
    (insert "denergy = " value))
  )


(defun pwcond-ecut2d ()
  (interactive)
  (let ((value (read-string "Value of ecut2d: ")))
    (insert "ecut2d = " value))
  )


(defun pwcond-energy0 ()
  (interactive)
  (let ((value (read-string "Value of energy0: ")))
    (insert "energy0 = " value))
  )


(defun pwcond-epsproj ()
  (interactive)
  (let ((value (read-string "Value of epsproj: ")))
    (insert "epsproj = " value))
  )


(defun pwcond-ewind ()
  (interactive)
  (let ((value (read-string "Value of ewind: ")))
    (insert "ewind = " value))
  )


(defun pwcond-fil_loc ()
  (interactive)
  (let ((value (read-string "Value of fil_loc: ")))
    (insert "fil_loc = '" value "'"))
  (backward-char 1)
  )


(defun pwcond-ikind ()
  (interactive)
  (let ((value (read-string "Value of ikind: ")))
    (insert "ikind = " value))
  )


(defun pwcond-iofspin ()
  (interactive)
  (let ((value (read-string "Value of iofspin: ")))
    (insert "iofspin = " value))
  )


(defun pwcond-last_e ()
  (interactive)
  (let ((value (read-string "Value of last_e: ")))
    (insert "last_e = " value))
  )


(defun pwcond-last_k ()
  (interactive)
  (let ((value (read-string "Value of last_k: ")))
    (insert "last_k = " value))
  )


(defun pwcond-llocal ()
  (interactive)
  (let ((value (read-string "Value of llocal: ")))
    (insert "llocal = " value))
  )


(defun pwcond-loop_ek ()
  (interactive)
  (let ((value (read-string "Value of loop_ek: ")))
    (insert "loop_ek = " value))
  )


(defun pwcond-lread_cond ()
  (interactive)
  (let ((value (read-string "Value of lread_cond: ")))
    (insert "lread_cond = " value))
  )


(defun pwcond-lread_loc ()
  (interactive)
  (let ((value (read-string "Value of lread_loc: ")))
    (insert "lread_loc = " value))
  )


(defun pwcond-lwrite_cond ()
  (interactive)
  (let ((value (read-string "Value of lwrite_cond: ")))
    (insert "lwrite_cond = " value))
  )


(defun pwcond-lwrite_loc ()
  (interactive)
  (let ((value (read-string "Value of lwrite_loc: ")))
    (insert "lwrite_loc = " value))
  )


(defun pwcond-max_seconds ()
  (interactive)
  (let ((value (read-string "Value of max_seconds: ")))
    (insert "max_seconds = " value))
  )


(defun pwcond-nenergy ()
  (interactive)
  (let ((value (read-string "Value of nenergy: ")))
    (insert "nenergy = " value))
  )


(defun pwcond-nz1 ()
  (interactive)
  (let ((value (read-string "Value of nz1: ")))
    (insert "nz1 = " value))
  )


(defun pwcond-orbj_fin ()
  (interactive)
  (let ((value (read-string "Value of orbj_fin: ")))
    (insert "orbj_fin = " value))
  )


(defun pwcond-orbj_in ()
  (interactive)
  (let ((value (read-string "Value of orbj_in: ")))
    (insert "orbj_in = " value))
  )


(defun pwcond-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun pwcond-prefixl ()
  (interactive)
  (let ((value (read-string "Value of prefixl: ")))
    (insert "prefixl = '" value "'"))
  (backward-char 1)
  )


(defun pwcond-prefixr ()
  (interactive)
  (let ((value (read-string "Value of prefixr: ")))
    (insert "prefixr = '" value "'"))
  (backward-char 1)
  )


(defun pwcond-prefixs ()
  (interactive)
  (let ((value (read-string "Value of prefixs: ")))
    (insert "prefixs = '" value "'"))
  (backward-char 1)
  )


(defun pwcond-prefixt ()
  (interactive)
  (let ((value (read-string "Value of prefixt: ")))
    (insert "prefixt = '" value "'"))
  (backward-char 1)
  )


(defun pwcond-recover ()
  (interactive)
  (let ((value (read-string "Value of recover: ")))
    (insert "recover = " value))
  )


(defun pwcond-save_file ()
  (interactive)
  (let ((value (read-file-name "Value of save_file: ")))
    (insert "save_file = '" value "'"))
  (backward-char 1)
  )


(defun pwcond-start_e ()
  (interactive)
  (let ((value (read-string "Value of start_e: ")))
    (insert "start_e = " value))
  )


(defun pwcond-start_k ()
  (interactive)
  (let ((value (read-string "Value of start_k: ")))
    (insert "start_k = " value))
  )


(defun pwcond-tk_plot ()
  (interactive)
  (let ((value (read-string "Value of tk_plot: ")))
    (insert "tk_plot = " value))
  )


(defun pwcond-tran_file ()
  (interactive)
  (let ((value (read-file-name "Value of tran_file: ")))
    (insert "tran_file = '" value "'"))
  (backward-char 1)
  )


(defun pwcond-tran_prefix ()
  (interactive)
  (let ((value (read-string "Value of tran_prefix: ")))
    (insert "tran_prefix = '" value "'"))
  (backward-char 1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; lanczos- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun lanczos-LR_CONTROL ()
  (interactive)
  (insert "&LR_CONTROL")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun lanczos-LR_INPUT ()
  (interactive)
  (insert "&LR_INPUT")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun lanczos-LR_POST ()
  (interactive)
  (insert "&LR_POST")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; lanczos- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun lanczos-beta_gamma_z_prefix ()
  (interactive)
  (let ((value (read-string "Value of beta_gamma_z_prefix: ")))
    (insert "beta_gamma_z_prefix = '" value "'"))
  (backward-char 1)
  )


(defun lanczos-charge_response ()
  (interactive)
  (let ((value (read-string "Value of charge_response: ")))
    (insert "charge_response = " value))
  )


(defun lanczos-d0psi_rs ()
  (interactive)
  (let ((value (read-string "Value of d0psi_rs: ")))
    (insert "d0psi_rs = " value))
  )


(defun lanczos-disk_io ()
  (interactive)
  (let ((value (read-string "Value of disk_io: ")))
    (insert "disk_io = '" value "'"))
  (backward-char 1)
  )


(defun lanczos-epsil ()
  (interactive)
  (let ((value (read-string "Value of epsil: ")))
    (insert "epsil = " value))
  )


(defun lanczos-ipol ()
  (interactive)
  (let ((value (read-string "Value of ipol: ")))
    (insert "ipol = " value))
  )


(defun lanczos-itermax ()
  (interactive)
  (let ((value (read-string "Value of itermax: ")))
    (insert "itermax = " value))
  )


(defun lanczos-lr_verbosity ()
  (interactive)
  (let ((value (read-string "Value of lr_verbosity: ")))
    (insert "lr_verbosity = " value))
  )


(defun lanczos-lrpa ()
  (interactive)
  (let ((value (read-string "Value of lrpa: ")))
    (insert "lrpa = " value))
  )


(defun lanczos-lshift_d0psi ()
  (interactive)
  (let ((value (read-string "Value of lshift_d0psi: ")))
    (insert "lshift_d0psi = " value))
  )


(defun lanczos-ltammd ()
  (interactive)
  (let ((value (read-string "Value of ltammd: ")))
    (insert "ltammd = " value))
  )


(defun lanczos-n_ipol ()
  (interactive)
  (let ((value (read-string "Value of n_ipol: ")))
    (insert "n_ipol = " value))
  )


(defun lanczos-no_hxc ()
  (interactive)
  (let ((value (read-string "Value of no_hxc: ")))
    (insert "no_hxc = " value))
  )


(defun lanczos-omeg ()
  (interactive)
  (let ((value (read-string "Value of omeg: ")))
    (insert "omeg = " value))
  )


(defun lanczos-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun lanczos-plot_type ()
  (interactive)
  (let ((value (read-string "Value of plot_type: ")))
    (insert "plot_type = " value))
  )


(defun lanczos-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun lanczos-pseudo_hermitian ()
  (interactive)
  (let ((value (read-string "Value of pseudo_hermitian: ")))
    (insert "pseudo_hermitian = " value))
  )


(defun lanczos-restart ()
  (interactive)
  (let ((value (read-string "Value of restart: ")))
    (insert "restart = " value))
  )


(defun lanczos-restart_step ()
  (interactive)
  (let ((value (read-string "Value of restart_step: ")))
    (insert "restart_step = " value))
  )


(defun lanczos-title ()
  (interactive)
  (let ((value (read-string "Value of title: ")))
    (insert "title = '" value "'"))
  (backward-char 1)
  )


(defun lanczos-w_t_npol ()
  (interactive)
  (let ((value (read-string "Value of w_T_npol: ")))
    (insert "w_T_npol = " value))
  )


(defun lanczos-wfcdir ()
  (interactive)
  (let ((value (read-directory-name "Value of wfcdir: ")))
    (insert "wfcdir = '" value "'"))
  (backward-char 1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; davidson- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun davidson-LR_DAV ()
  (interactive)
  (insert "&LR_DAV")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun davidson-LR_INPUT ()
  (interactive)
  (insert "&LR_INPUT")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; davidson- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun davidson-broadening ()
  (interactive)
  (let ((value (read-string "Value of broadening: ")))
    (insert "broadening = " value))
  )


(defun davidson-d0psi_rs ()
  (interactive)
  (let ((value (read-string "Value of d0psi_rs: ")))
    (insert "d0psi_rs = " value))
  )


(defun davidson-disk_io ()
  (interactive)
  (let ((value (read-string "Value of disk_io: ")))
    (insert "disk_io = '" value "'"))
  (backward-char 1)
  )


(defun davidson-finish ()
  (interactive)
  (let ((value (read-string "Value of finish: ")))
    (insert "finish = " value))
  )


(defun davidson-if_dft_spectrum ()
  (interactive)
  (let ((value (read-string "Value of if_dft_spectrum: ")))
    (insert "if_dft_spectrum = " value))
  )


(defun davidson-if_random_init ()
  (interactive)
  (let ((value (read-string "Value of if_random_init: ")))
    (insert "if_random_init = " value))
  )


(defun davidson-lplot_drho ()
  (interactive)
  (let ((value (read-string "Value of lplot_drho: ")))
    (insert "lplot_drho = " value))
  )


(defun davidson-lr_verbosity ()
  (interactive)
  (let ((value (read-string "Value of lr_verbosity: ")))
    (insert "lr_verbosity = " value))
  )


(defun davidson-lshift_d0psi ()
  (interactive)
  (let ((value (read-string "Value of lshift_d0psi: ")))
    (insert "lshift_d0psi = " value))
  )


(defun davidson-ltammd ()
  (interactive)
  (let ((value (read-string "Value of ltammd: ")))
    (insert "ltammd = " value))
  )


(defun davidson-max_iter ()
  (interactive)
  (let ((value (read-string "Value of max_iter: ")))
    (insert "max_iter = " value))
  )


(defun davidson-max_seconds ()
  (interactive)
  (let ((value (read-string "Value of max_seconds: ")))
    (insert "max_seconds = " value))
  )


(defun davidson-no_hxc ()
  (interactive)
  (let ((value (read-string "Value of no_hxc: ")))
    (insert "no_hxc = " value))
  )


(defun davidson-num_basis_max ()
  (interactive)
  (let ((value (read-string "Value of num_basis_max: ")))
    (insert "num_basis_max = " value))
  )


(defun davidson-num_eign ()
  (interactive)
  (let ((value (read-string "Value of num_eign: ")))
    (insert "num_eign = " value))
  )


(defun davidson-num_init ()
  (interactive)
  (let ((value (read-string "Value of num_init: ")))
    (insert "num_init = " value))
  )


(defun davidson-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun davidson-p_nbnd_occ ()
  (interactive)
  (let ((value (read-string "Value of p_nbnd_occ: ")))
    (insert "p_nbnd_occ = " value))
  )


(defun davidson-p_nbnd_virt ()
  (interactive)
  (let ((value (read-string "Value of p_nbnd_virt: ")))
    (insert "p_nbnd_virt = " value))
  )


(defun davidson-poor_of_ram ()
  (interactive)
  (let ((value (read-string "Value of poor_of_ram: ")))
    (insert "poor_of_ram = " value))
  )


(defun davidson-poor_of_ram2 ()
  (interactive)
  (let ((value (read-string "Value of poor_of_ram2: ")))
    (insert "poor_of_ram2 = " value))
  )


(defun davidson-precondition ()
  (interactive)
  (let ((value (read-string "Value of precondition: ")))
    (insert "precondition = " value))
  )


(defun davidson-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun davidson-pseudo_hermitian ()
  (interactive)
  (let ((value (read-string "Value of pseudo_hermitian: ")))
    (insert "pseudo_hermitian = " value))
  )


(defun davidson-reference ()
  (interactive)
  (let ((value (read-string "Value of reference: ")))
    (insert "reference = " value))
  )


(defun davidson-residue_conv_thr ()
  (interactive)
  (let ((value (read-string "Value of residue_conv_thr: ")))
    (insert "residue_conv_thr = " value))
  )


(defun davidson-restart ()
  (interactive)
  (let ((value (read-string "Value of restart: ")))
    (insert "restart = " value))
  )


(defun davidson-single_pole ()
  (interactive)
  (let ((value (read-string "Value of single_pole: ")))
    (insert "single_pole = " value))
  )


(defun davidson-start ()
  (interactive)
  (let ((value (read-string "Value of start: ")))
    (insert "start = " value))
  )


(defun davidson-step ()
  (interactive)
  (let ((value (read-string "Value of step: ")))
    (insert "step = " value))
  )


(defun davidson-wfcdir ()
  (interactive)
  (let ((value (read-directory-name "Value of wfcdir: ")))
    (insert "wfcdir = '" value "'"))
  (backward-char 1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; spectrum- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun spectrum-LR_INPUT ()
  (interactive)
  (insert "&LR_INPUT")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; spectrum- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun spectrum-eels ()
  (interactive)
  (let ((value (read-string "Value of eels: ")))
    (insert "eels = " value))
  )


(defun spectrum-eign_file ()
  (interactive)
  (let ((value (read-file-name "Value of eign_file: ")))
    (insert "eign_file = '" value "'"))
  (backward-char 1)
  )


(defun spectrum-end ()
  (interactive)
  (let ((value (read-string "Value of end: ")))
    (insert "end = " value))
  )


(defun spectrum-epsil ()
  (interactive)
  (let ((value (read-string "Value of epsil: ")))
    (insert "epsil = " value))
  )


(defun spectrum-extrapolation ()
  (interactive)
  (let ((value (read-string "Value of extrapolation: ")))
    (insert "extrapolation = '" value "'"))
  (backward-char 1)
  )


(defun spectrum-increment ()
  (interactive)
  (let ((value (read-string "Value of increment: ")))
    (insert "increment = " value))
  )


(defun spectrum-ipol ()
  (interactive)
  (let ((value (read-string "Value of ipol: ")))
    (insert "ipol = " value))
  )


(defun spectrum-itermax ()
  (interactive)
  (let ((value (read-string "Value of itermax: ")))
    (insert "itermax = " value))
  )


(defun spectrum-itermax0 ()
  (interactive)
  (let ((value (read-string "Value of itermax0: ")))
    (insert "itermax0 = " value))
  )


(defun spectrum-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun spectrum-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun spectrum-start ()
  (interactive)
  (let ((value (read-string "Value of start: ")))
    (insert "start = " value))
  )


(defun spectrum-td ()
  (interactive)
  (let ((value (read-string "Value of td: ")))
    (insert "td = '" value "'"))
  (backward-char 1)
  )


(defun spectrum-units ()
  (interactive)
  (let ((value (read-string "Value of units: ")))
    (insert "units = " value))
  )


(defun spectrum-verbosity ()
  (interactive)
  (let ((value (read-string "Value of verbosity: ")))
    (insert "verbosity = " value))
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; eels- namelists functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun eels-LR_CONTROL ()
  (interactive)
  (insert "&LR_CONTROL")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )


(defun eels-LR_INPUT ()
  (interactive)
  (insert "&LR_INPUT")
  (newline 2)
  (insert "/")
  (forward-line -1)
  )



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; eels- namelist's variables functions ...
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun eels-approximation ()
  (interactive)
  (let ((value (read-string "Value of approximation: ")))
    (insert "approximation = '" value "'"))
  (backward-char 1)
  )


(defun eels-disk_io ()
  (interactive)
  (let ((value (read-string "Value of disk_io: ")))
    (insert "disk_io = '" value "'"))
  (backward-char 1)
  )


(defun eels-itermax ()
  (interactive)
  (let ((value (read-string "Value of itermax: ")))
    (insert "itermax = " value))
  )


(defun eels-lr_verbosity ()
  (interactive)
  (let ((value (read-string "Value of lr_verbosity: ")))
    (insert "lr_verbosity = " value))
  )


(defun eels-outdir ()
  (interactive)
  (let ((value (read-directory-name "Value of outdir: ")))
    (insert "outdir = '" value "'"))
  (backward-char 1)
  )


(defun eels-prefix ()
  (interactive)
  (let ((value (read-string "Value of prefix: ")))
    (insert "prefix = '" value "'"))
  (backward-char 1)
  )


(defun eels-pseudo_hermitian ()
  (interactive)
  (let ((value (read-string "Value of pseudo_hermitian: ")))
    (insert "pseudo_hermitian = " value))
  )


(defun eels-qi ()
  (interactive)
  (let ((value (read-string "Value of qi: ")))
    (insert "qi = " value))
  )


(defun eels-restart ()
  (interactive)
  (let ((value (read-string "Value of restart: ")))
    (insert "restart = " value))
  )


(defun eels-restart_step ()
  (interactive)
  (let ((value (read-string "Value of restart_step: ")))
    (insert "restart_step = " value))
  )





(provide 'qe-funcs)



