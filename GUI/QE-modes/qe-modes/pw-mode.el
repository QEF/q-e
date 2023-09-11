;; pw-mode.el
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
;; (http://ergoemacs.org/emacs/elisp_syntax_coloring.html). Sebastijan
;; Peljhan is acknowledged for his work on `xsf-mode' that inspired
;; the idea of writing the qe-modes. Last but not the least,
;; Hongyi Zhao contributed the ido-completion-read snippet of
;; code for selecting the values for the card's flags.


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

;; This is the `pw-mode', a major mode for composing the Quantum ESPRESSO
;; QE-pw.x input files. For the installation and usage, see the
;; user_guide.pdf in the Doc/ subdirectory of the original package
;; (quick installation instructions are also available in the README
;; file of the original package).

;;; Code:

(require 'font-lock)
(require 'regexp-opt)

(defvar pw-mode-hook nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; basic variables

;; pw's supercards (if any)
(defvar pw-open-supercards   (list ))
(defvar pw-closed-supercards (list ))
  
;; pw's namelists
(defvar pw-namelists (list "&CELL" "&CONTROL" "&ELECTRONS" "&FCP" "&IONS" "&RISM" "&SYSTEM" ))
(defvar qe-end-namelist (list "&END" "/"))

;; pw's variables
(defvar pw-vars (list "A" "B" "C" "Hubbard_alpha" "Hubbard_beta" "Hubbard_occ" "ace" "adaptive_thr" "angle1" "angle2" "assume_isolated" "bfgs_ndim" "block" "block_1" "block_2" "block_height" "calculation" "cell_dofree" "cell_dynamics" "cell_factor" "celldm" "closure" "constrained_magnetization" "conv_thr" "conv_thr_init" "conv_thr_multi" "cosAB" "cosAC" "cosBC" "degauss" "degauss_cond" "delta_t" "dftd3_threebody" "dftd3_version" "diago_cg_maxiter" "diago_david_ndim" "diago_full_acc" "diago_gs_nblock" "diago_ppcg_maxiter" "diago_rmm_conv" "diago_rmm_ndim" "diago_thr_init" "diagonalization" "dipfield" "disk_io" "dmft" "dmft_prefix" "dt" "eamp" "ecfixed" "ecutfock" "ecutrho" "ecutsolv" "ecutvcut" "ecutwfc" "edir" "efield" "efield_cart" "efield_phase" "electron_maxstep" "emaxpos" "ensemble_energies" "eopreg" "esm_bc" "esm_efield" "esm_nfit" "esm_w" "etot_conv_thr" "exx_fraction" "exx_maxstep" "exxdiv_treatment" "fcp_conv_thr" "fcp_delta_t" "fcp_dynamics" "fcp_mass" "fcp_mu" "fcp_ndiis" "fcp_nraise" "fcp_temperature" "fcp_tempw" "fcp_tolp" "fcp_velocity" "fire_alpha_init" "fire_dtmax" "fire_f_dec" "fire_f_inc" "fire_falpha" "fire_nmin" "fixed_magnetization" "forc_conv_thr" "force_symmorphic" "freeze_all_atoms" "gate" "gcscf_beta" "gcscf_conv_thr" "gcscf_mu" "gdir" "ibrav" "input_dft" "ion_dynamics" "ion_positions" "ion_temperature" "ion_velocities" "iprint" "lambda" "laue_both_hands" "laue_buffer_left" "laue_buffer_right" "laue_expand_left" "laue_expand_right" "laue_nfit" "laue_starting_left" "laue_starting_right" "laue_wall" "laue_wall_epsilon" "laue_wall_lj6" "laue_wall_rho" "laue_wall_sigma" "laue_wall_z" "lberry" "lelfield" "lfcp" "lforcet" "lgcscf" "lkpoint_dir" "localization_thr" "london" "london_c6" "london_rcut" "london_rvdw" "london_s6" "lorbm" "lspinorb" "max_seconds" "mdiis1d_size" "mdiis1d_step" "mdiis3d_size" "mdiis3d_step" "mixing_beta" "mixing_fixed_ns" "mixing_mode" "mixing_ndim" "nat" "nberrycyc" "nbnd" "nbnd_cond" "nelec_cond" "nextffield" "no_t_rev" "noinv" "noncolin" "nosym" "nosym_evc" "nppstr" "nqx1" "nqx2" "nqx3" "nr1" "nr1s" "nr2" "nr2s" "nr3" "nr3s" "nraise" "nsolv" "nspin" "nstep" "ntyp" "occupations" "one_atom_occupations" "origin_choice" "outdir" "pol_type" "pot_extrapolation" "prefix" "press" "press_conv_thr" "pseudo_dir" "q2sigma" "qcutz" "real_space" "refold_pos" "relaxz" "remove_rigid_rot" "report" "restart_mode" "rhombohedral" "rism1d_bond_width" "rism1d_conv_thr" "rism1d_dielectric" "rism1d_maxstep" "rism1d_molesize" "rism1d_nproc" "rism3d_conv_level" "rism3d_conv_thr" "rism3d_maxstep" "rism3d_planar_average" "scf_must_converge" "sci_cb" "sci_vb" "screening_parameter" "sic_energy" "sic_gamma" "smear1d" "smear3d" "smearing" "solute_epsilon" "solute_lj" "solute_sigma" "space_group" "starting1d" "starting3d" "starting_charge" "starting_magnetization" "starting_ns_eigenvalue" "starting_spin_angle" "startingpot" "startingwfc" "tefield" "tempv" "tempw" "title" "tolp" "tot_charge" "tot_magnetization" "tprnfor" "tqr" "trism" "trust_radius_ini" "trust_radius_max" "trust_radius_min" "ts_vdw_econv_thr" "ts_vdw_isolated" "tstress" "twochem" "uniqueb" "upscale" "use_all_frac" "vdw_corr" "verbosity" "w_1" "w_2" "wf_collect" "wfc_extrapolation" "wfcdir" "wmass" "x_gamma_extrapolation" "xdm" "xdm_a1" "xdm_a2" "zgate" ))

;; pw's cards & keywords
(defvar pw-cards (list "ADDITIONAL_K_POINTS" "ATOMIC_FORCES" "ATOMIC_POSITIONS" "ATOMIC_SPECIES" "ATOMIC_VELOCITIES" "CELL_PARAMETERS" "CONSTRAINTS" "HUBBARD" "J0" "K_POINTS" "OCCUPATIONS" "SOLVENTS" "U" "V" ))

;; pw's flags
(defvar pw-flags (list "1/cell" "a.u" "alat" "angstrom" "atomic" "automatic" "bohr" "crystal" "crystal_b" "crystal_c" "crystal_sg" "g/cm^3" "gamma" "mol/L" "norm-atomic" "ortho-atomic" "pseudo" "tpiba" "tpiba_b" "tpiba_c" "wf" ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; derived variables
  
(defvar pw-open-supercards-regexp   (regexp-opt pw-open-supercards   'symbols)) ; may not exists
(defvar pw-closed-supercards-regexp (regexp-opt pw-closed-supercards 'symbols)) ; may not exists

(defvar pw-cards-regexp (regexp-opt
			    (append pw-cards pw-open-supercards) 'symbols))
(defvar pw-flags-regexp (regexp-opt pw-flags 'symbols))

(defvar pw-namelist-face (cons (regexp-opt (append pw-namelists qe-end-namelist) 'symbols) font-lock-function-name-face))
(defvar pw-variable-face (cons (regexp-opt pw-vars 'symbols) font-lock-variable-name-face))

;; logical values as constants
(defvar qe-logic-face (cons (regexp-opt (list ".t." ".true." ".f." ".false.")) font-lock-constant-face))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; regexp for indentation
(defvar pw-decr-indent-fold-t-re (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
(defvar pw-decr-indent-re        (concat "^[ \t]*" (regexp-opt
						       (append pw-cards pw-open-supercards pw-closed-supercards) t)))
;;
(defvar pw-deindent-fold-t-re    (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
;;
(defvar pw-indent-fold-t-re      (concat "^[ \t]*" (regexp-opt pw-namelists t)))
(defvar pw-indent-re             (concat "^[ \t]*" (regexp-opt pw-cards     t)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; supercards, cards and flags are case sensitive -- here are the corresponding matchers

(defun pw-closed-supercards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward pw-closed-supercards-regexp limit 'no-error)))

(defun pw-cards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward pw-cards-regexp limit 'no-error)))

(defun pw-flags-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward pw-flags-regexp limit 'no-error)))


(font-lock-add-keywords 'pw-mode (list
				     pw-namelist-face 
				     pw-variable-face
				     qe-logic-face
				     '("," . font-lock-builtin-face)
				     '("(" . font-lock-builtin-face)
				     '(")" . font-lock-builtin-face)
				     ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; register the keywords

(font-lock-add-keywords 'pw-mode '(
				      (pw-closed-supercards-matcher 1 font-lock-preprocessor-face t)
				      (pw-cards-matcher 1 font-lock-keyword-face t)
				      (pw-flags-matcher 1 font-lock-type-face    t)
				      ))

;;(defvar pw-keywords '(pw-namelist-face pw-variable-face))
(defvar pw-keywords '(((list "") . font-lock-constant-face)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the sytnax of strings

(defvar pw-mode-syntax-table
  (let ((table (make-syntax-table)))
    (modify-syntax-entry ?\' "\"'"  table)
    (modify-syntax-entry ?\" "\"\"" table)
    table)
  "Syntax table in use in `pw-mode' buffers.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; code for auto-indenting

(defvar qe-indent 3)
(defun pw-indent-line ()
  "Indent current line according to pw input syntax."
  (interactive)
  (beginning-of-line)
  (if (bobp)
      (indent-line-to 0)		   ; First line indented to column 0
    (let ((not-indented t) cur-indent)
      (if (or (looking-at pw-decr-indent-fold-t-re)
	      (let ((case-fold-search nil)) (looking-at pw-decr-indent-re))) ; If the line we are looking at is the end of a block, then decrease the indentation
	  (progn
	    (save-excursion
	      (forward-line -1)
	      (setq cur-indent (- (current-indentation) qe-indent)))
	    (if (< cur-indent 0) ; We can't indent past the left margin
		(setq cur-indent 0)))
	(save-excursion
	  (while  not-indented ; Iterate backwards until we find an indentation hint
	    (forward-line -1)
	    (if (looking-at pw-deindent-fold-t-re) ; This hint indicates that we need to indent at the level of the "/" token
		(progn
		  (setq cur-indent (current-indentation))
		  (setq not-indented nil))
	      (if (or (looking-at pw-indent-fold-t-re)
		      (let ((case-fold-search nil)) (looking-at pw-indent-re))) ; This hint indicates that we need to indent an extra level
		  (progn
		    (setq cur-indent (+ (current-indentation) qe-indent)) ; Do the actual indenting
		    (setq not-indented nil))
		(if (bobp)
		    (setq not-indented nil)))))))
      (if cur-indent
	  (indent-line-to cur-indent)
	(indent-line-to 0))))) ; If we didn't see an indentation hint, then allow no indentation


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the pw-mode as derived-mode

(define-derived-mode pw-mode prog-mode 
  "QE-pw.x"  
  "Major mode for editing Qunatum-ESPRESSO input files (QE-pw.x mode)"  
  (setq font-lock-defaults '((pw-keywords) nil t))
  (set (make-local-variable 'indent-line-function) 'pw-indent-line)
  
  ;; define the syntax of comments
  (setq comment-start "!")
  (setq comment-end "")
  (modify-syntax-entry ?!  "< b" pw-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" pw-mode-syntax-table)
  (modify-syntax-entry ?=  " " pw-mode-syntax-table) ;; treat "=" non symbol constituent
  ;; end
  )

;; free memory

(setq pw-namelists nil)
(setq pw-vars nil)
(setq pw-cards nil)
(setq pw-flags nil)
(setq pw-open-supercards   nil)
(setq pw-closed-supercards nil)


(require 'qe-funcs)

(provide 'pw-mode)

