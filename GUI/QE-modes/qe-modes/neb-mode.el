;; neb-mode.el
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

;; This is the `neb-mode', a major mode for composing the Quantum ESPRESSO
;; QE-neb.x input files. For the installation and usage, see the
;; user_guide.pdf in the Doc/ subdirectory of the original package
;; (quick installation instructions are also available in the README
;; file of the original package).

;;; Code:

(require 'font-lock)
(require 'regexp-opt)

(defvar neb-mode-hook nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; basic variables

;; neb's supercards (if any)
(defvar neb-open-supercards   (list "FIRST_IMAGE" "INTERMEDIATE_IMAGE" "LAST_IMAGE" ))
(defvar neb-closed-supercards (list "BEGIN" "BEGIN_ENGINE_INPUT" "BEGIN_PATH_INPUT" "BEGIN_POSITIONS" "END" "END_ENGINE_INPUT" "END_PATH_INPUT" "END_POSITIONS" ))
  
;; neb's namelists
(defvar neb-namelists (list "&CELL" "&CONTROL" "&ELECTRONS" "&IONS" "&PATH" "&SYSTEM" ))
(defvar qe-end-namelist (list "&END" "/"))

;; neb's variables
(defvar neb-vars (list "A" "adaptive_thr" "angle1" "angle2" "assume_isolated" "B" "bfgs_ndim" "block" "block_1" "block_2" "block_height" "C" "calculation" "cell_dofree" "cell_dynamics" "cell_factor" "celldm" "CI_scheme" "constrained_magnetization" "conv_thr" "conv_thr_init" "conv_thr_multi" "cosAB" "cosAC" "cosBC" "degauss" "delta_t" "dftd3_threebody" "dftd3_version" "diago_cg_maxiter" "diago_david_ndim" "diago_full_acc" "diago_thr_init" "diagonalization" "dipfield" "disk_io" "ds" "dt" "eamp" "ecfixed" "ecutfock" "ecutrho" "ecutvcut" "ecutwfc" "edir" "efield" "efield_cart" "efield_phase" "electron_maxstep" "emaxpos" "eopreg" "esm_bc" "esm_efield" "esm_nfit" "esm_w" "etot_conv_thr" "exx_fraction" "exxdiv_treatment" "fcp_mu" "fcp_tot_charge_first" "fcp_tot_charge_last" "first_last_opt" "fixed_magnetization" "forc_conv_thr" "force_symmorphic" "gate" "gdir" "Hubbard_alpha" "Hubbard_beta" "Hubbard_J" "Hubbard_J0" "Hubbard_U" "ibrav" "input_dft" "ion_dynamics" "ion_positions" "ion_temperature" "iprint" "k_max" "k_min" "lambda" "lberry" "lda_plus_u" "lda_plus_u_kind" "lelfield" "lfcpopt" "lforcet" "lkpoint_dir" "localization_thr" "london" "london_c6" "london_rcut" "london_rvdw" "london_s6" "lorbm" "lspinorb" "max_seconds" "minimum_image" "mixing_beta" "mixing_fixed_ns" "mixing_mode" "mixing_ndim" "nat" "nberrycyc" "nbnd" "no_t_rev" "noinv" "noncolin" "nosym" "nosym_evc" "nppstr" "nqx1" "nqx2" "nqx3" "nr1" "nr1s" "nr2" "nr2s" "nr3" "nr3s" "nraise" "nspin" "nstep" "nstep_path" "ntyp" "num_of_images" "occupations" "one_atom_occupations" "opt_scheme" "origin_choice" "ortho_para" "outdir" "path_thr" "pot_extrapolation" "prefix" "press" "press_conv_thr" "pseudo_dir" "q2sigma" "qcutz" "real_space" "refold_pos" "relaxz" "remove_rigid_rot" "report" "restart_mode" "rhombohedral" "scf_must_converge" "screening_parameter" "smearing" "space_group" "starting_charge" "starting_magnetization" "starting_ns_eigenvalue" "starting_spin_angle" "startingpot" "startingwfc" "string_method" "tefield" "temp_req" "tempw" "title" "tolp" "tot_charge" "tot_magnetization" "tprnfor" "tqr" "trust_radius_ini" "trust_radius_max" "trust_radius_min" "ts_vdw_econv_thr" "ts_vdw_isolated" "tstress" "U_projection_type" "uniqueb" "upscale" "use_all_frac" "use_freezing" "use_masses" "vdw_corr" "verbosity" "w_1" "w_2" "wf_collect" "wfc_extrapolation" "wfcdir" "wmass" "x_gamma_extrapolation" "xdm" "xdm_a1" "xdm_a2" "zgate" ))

;; neb's cards & keywords
(defvar neb-cards (list "ATOMIC_FORCES" "ATOMIC_POSITIONS" "ATOMIC_SPECIES" "CELL_PARAMETERS" "CLIMBING_IMAGES" "CONSTRAINTS" "K_POINTS" "OCCUPATIONS" ))

;; neb's flags
(defvar neb-flags (list "alat" "angstrom" "automatic" "bohr" "crystal" "crystal_b" "crystal_c" "crystal_sg" "gamma" "tpiba" "tpiba_b" "tpiba_c" ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; derived variables
  
(defvar neb-open-supercards-regexp   (regexp-opt neb-open-supercards   'symbols)) ; may not exists
(defvar neb-closed-supercards-regexp (regexp-opt neb-closed-supercards 'symbols)) ; may not exists

(defvar neb-cards-regexp (regexp-opt
			    (append neb-cards neb-open-supercards) 'symbols))
(defvar neb-flags-regexp (regexp-opt neb-flags 'symbols))

(defvar neb-namelist-face (cons (regexp-opt (append neb-namelists qe-end-namelist) 'symbols) font-lock-function-name-face))
(defvar neb-variable-face (cons (regexp-opt neb-vars 'symbols) font-lock-variable-name-face))

;; logical values as constants
(defvar qe-logic-face (cons (regexp-opt (list ".t." ".true." ".f." ".false.")) font-lock-constant-face))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; regexp for indentation
(defvar neb-decr-indent-fold-t-re (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
(defvar neb-decr-indent-re        (concat "^[ \t]*" (regexp-opt
						       (append neb-cards neb-open-supercards neb-closed-supercards) t)))
;;
(defvar neb-deindent-fold-t-re    (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
;;
(defvar neb-indent-fold-t-re      (concat "^[ \t]*" (regexp-opt neb-namelists t)))
(defvar neb-indent-re             (concat "^[ \t]*" (regexp-opt neb-cards     t)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; supercards, cards and flags are case sensitive -- here are the corresponding matchers

(defun neb-closed-supercards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward neb-closed-supercards-regexp limit 'no-error)))

(defun neb-cards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward neb-cards-regexp limit 'no-error)))

(defun neb-flags-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward neb-flags-regexp limit 'no-error)))


(font-lock-add-keywords 'neb-mode (list
				     neb-namelist-face 
				     neb-variable-face
				     qe-logic-face
				     '("," . font-lock-builtin-face)
				     '("(" . font-lock-builtin-face)
				     '(")" . font-lock-builtin-face)
				     ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; register the keywords

(font-lock-add-keywords 'neb-mode '(
				      (neb-closed-supercards-matcher 1 font-lock-preprocessor-face t)
				      (neb-cards-matcher 1 font-lock-keyword-face t)
				      (neb-flags-matcher 1 font-lock-type-face    t)
				      ))

;;(defvar neb-keywords '(neb-namelist-face neb-variable-face))
(defvar neb-keywords '(((list "") . font-lock-constant-face)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the sytnax of strings

(defvar neb-mode-syntax-table
  (let ((table (make-syntax-table)))
    (modify-syntax-entry ?\' "\"'"  table)
    (modify-syntax-entry ?\" "\"\"" table)
    table)
  "Syntax table in use in `neb-mode' buffers.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; code for auto-indenting

(defvar qe-indent 3)
(defun neb-indent-line ()
  "Indent current line according to neb input syntax."
  (interactive)
  (beginning-of-line)
  (if (bobp)
      (indent-line-to 0)		   ; First line indented to column 0
    (let ((not-indented t) cur-indent)
      (if (or (looking-at neb-decr-indent-fold-t-re)
	      (let ((case-fold-search nil)) (looking-at neb-decr-indent-re))) ; If the line we are looking at is the end of a block, then decrease the indentation
	  (progn
	    (save-excursion
	      (forward-line -1)
	      (setq cur-indent (- (current-indentation) qe-indent)))
	    (if (< cur-indent 0) ; We can't indent past the left margin
		(setq cur-indent 0)))
	(save-excursion
	  (while  not-indented ; Iterate backwards until we find an indentation hint
	    (forward-line -1)
	    (if (looking-at neb-deindent-fold-t-re) ; This hint indicates that we need to indent at the level of the "/" token
		(progn
		  (setq cur-indent (current-indentation))
		  (setq not-indented nil))
	      (if (or (looking-at neb-indent-fold-t-re)
		      (let ((case-fold-search nil)) (looking-at neb-indent-re))) ; This hint indicates that we need to indent an extra level
		  (progn
		    (setq cur-indent (+ (current-indentation) qe-indent)) ; Do the actual indenting
		    (setq not-indented nil))
		(if (bobp)
		    (setq not-indented nil)))))))
      (if cur-indent
	  (indent-line-to cur-indent)
	(indent-line-to 0))))) ; If we didn't see an indentation hint, then allow no indentation


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the neb-mode as derived-mode

(define-derived-mode neb-mode prog-mode 
  "QE-neb.x"  
  "Major mode for editing Qunatum-ESPRESSO input files (QE-neb.x mode)"  
  (setq font-lock-defaults '((neb-keywords) nil t))
  (set (make-local-variable 'indent-line-function) 'neb-indent-line)
  
  ;; define the syntax of comments
  (setq comment-start "!")
  (setq comment-end "")
  (modify-syntax-entry ?!  "< b" neb-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" neb-mode-syntax-table)
  (modify-syntax-entry ?=  " " neb-mode-syntax-table) ;; treat "=" non symbol constituent
  ;; end
  )

;; free memory

(setq neb-namelists nil)
(setq neb-vars nil)
(setq neb-cards nil)
(setq neb-flags nil)
(setq neb-open-supercards   nil)
(setq neb-closed-supercards nil)


(require 'qe-funcs)

(provide 'neb-mode)

