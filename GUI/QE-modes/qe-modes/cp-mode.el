;; cp-mode.el
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

;; This is the `cp-mode', a major mode for composing the Quantum ESPRESSO
;; QE-cp.x input files. For the installation and usage, see the
;; user_guide.pdf in the Doc/ subdirectory of the original package
;; (quick installation instructions are also available in the README
;; file of the original package).

;;; Code:

(require 'font-lock)
(require 'regexp-opt)

(defvar cp-mode-hook nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; basic variables

;; cp's supercards (if any)
(defvar cp-open-supercards   (list ))
(defvar cp-closed-supercards (list ))
  
;; cp's namelists
(defvar cp-namelists (list "&CELL" "&CONTROL" "&ELECTRONS" "&IONS" "&PRESS_AI" "&SYSTEM" "&WANNIER" ))
(defvar qe-end-namelist (list "&END" "/"))

;; cp's variables
(defvar cp-vars (list "A" "abisur" "abivol" "adapt" "ampre" "amprp" "assume_isolated" "B" "C" "calculation" "calwf" "cell_damping" "cell_dofree" "cell_dynamics" "cell_factor" "cell_parameters" "cell_temperature" "cell_velocities" "celldm" "conv_thr" "cosAB" "cosAC" "cosBC" "degauss" "disk_io" "dt" "dthr" "ecfixed" "ecutrho" "ecutwfc" "efield" "efx0" "efx1" "efy0" "efy1" "efz0" "efz1" "ekin_conv_thr" "ekincw" "electron_damping" "electron_dynamics" "electron_maxstep" "electron_temperature" "electron_velocities" "emass" "emass_cutoff" "epol" "etot_conv_thr" "exx_dis_cutoff" "exx_fraction" "exx_me_rcut_pair" "exx_me_rcut_self" "exx_neigh" "exx_poisson_eps" "exx_ps_rcut_pair" "exx_ps_rcut_self" "fnhscl" "fnosee" "fnoseh" "fnosep" "forc_conv_thr" "grease" "greash" "greasp" "Hubbard_U" "ibrav" "iesr" "input_dft" "ion_damping" "ion_dynamics" "ion_nstepe" "ion_positions" "ion_radius" "ion_temperature" "ion_velocities" "iprint" "isave" "lambda_cold" "lda_plus_u" "london_rcut" "london_s6" "max_seconds" "maxiter" "maxwfdt" "memory" "n_inner" "nat" "nbnd" "ndega" "ndr" "ndw" "nhgrp" "nhpcl" "nhptyp" "ninter_cold_restart" "nit" "niter_cg_restart" "nr1" "nr1b" "nr1s" "nr2" "nr2b" "nr2s" "nr3" "nr3b" "nr3s" "nsd" "nspin" "nstep" "nsteps" "ntyp" "nwf" "occupations" "ortho_eps" "ortho_max" "ortho_para" "orthogonalization" "outdir" "P_ext" "P_fin" "P_in" "passop" "prefix" "press" "pseudo_dir" "pvar" "q2sigma" "qcutz" "remove_rigid_rot" "restart_mode" "rho_thr" "saverho" "smearing" "startingwfc" "Surf_t" "sw_len" "tabps" "tcg" "tefield" "temph" "tempw" "title" "tolp" "tolw" "tot_charge" "tot_magnetization" "tprnfor" "tranp" "ts_vdw" "ts_vdw_econv_thr" "ts_vdw_isolated" "tstress" "vdw_corr" "verbosity" "wf_efield" "wf_friction" "wf_q" "wf_switch" "wfdt" "wffort" "wfsd" "wmass" "writev" ))

;; cp's cards & keywords
(defvar cp-cards (list "ATOMIC_FORCES" "ATOMIC_POSITIONS" "ATOMIC_SPECIES" "ATOMIC_VELOCITIES" "AUTOPILOT" "CELL_PARAMETERS" "CONSTRAINTS" "OCCUPATIONS" "PLOT_WANNIER" "REF_CELL_PARAMETERS" ))

;; cp's flags
(defvar cp-flags (list "a.u" "alat" "angstrom" "bohr" "crystal" ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; derived variables
  
(defvar cp-open-supercards-regexp   (regexp-opt cp-open-supercards   'symbols)) ; may not exists
(defvar cp-closed-supercards-regexp (regexp-opt cp-closed-supercards 'symbols)) ; may not exists

(defvar cp-cards-regexp (regexp-opt
			    (append cp-cards cp-open-supercards) 'symbols))
(defvar cp-flags-regexp (regexp-opt cp-flags 'symbols))

(defvar cp-namelist-face (cons (regexp-opt (append cp-namelists qe-end-namelist) 'symbols) font-lock-function-name-face))
(defvar cp-variable-face (cons (regexp-opt cp-vars 'symbols) font-lock-variable-name-face))

;; logical values as constants
(defvar qe-logic-face (cons (regexp-opt (list ".t." ".true." ".f." ".false.")) font-lock-constant-face))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; regexp for indentation
(defvar cp-decr-indent-fold-t-re (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
(defvar cp-decr-indent-re        (concat "^[ \t]*" (regexp-opt
						       (append cp-cards cp-open-supercards cp-closed-supercards) t)))
;;
(defvar cp-deindent-fold-t-re    (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
;;
(defvar cp-indent-fold-t-re      (concat "^[ \t]*" (regexp-opt cp-namelists t)))
(defvar cp-indent-re             (concat "^[ \t]*" (regexp-opt cp-cards     t)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; supercards, cards and flags are case sensitive -- here are the corresponding matchers

(defun cp-closed-supercards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward cp-closed-supercards-regexp limit 'no-error)))

(defun cp-cards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward cp-cards-regexp limit 'no-error)))

(defun cp-flags-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward cp-flags-regexp limit 'no-error)))


(font-lock-add-keywords 'cp-mode (list
				     cp-namelist-face 
				     cp-variable-face
				     qe-logic-face
				     '("," . font-lock-builtin-face)
				     '("(" . font-lock-builtin-face)
				     '(")" . font-lock-builtin-face)
				     ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; register the keywords

(font-lock-add-keywords 'cp-mode '(
				      (cp-closed-supercards-matcher 1 font-lock-preprocessor-face t)
				      (cp-cards-matcher 1 font-lock-keyword-face t)
				      (cp-flags-matcher 1 font-lock-type-face    t)
				      ))

;;(defvar cp-keywords '(cp-namelist-face cp-variable-face))
(defvar cp-keywords '(((list "") . font-lock-constant-face)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the sytnax of strings

(defvar cp-mode-syntax-table
  (let ((table (make-syntax-table)))
    (modify-syntax-entry ?\' "\"'"  table)
    (modify-syntax-entry ?\" "\"\"" table)
    table)
  "Syntax table in use in `cp-mode' buffers.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; code for auto-indenting

(defvar qe-indent 3)
(defun cp-indent-line ()
  "Indent current line according to cp input syntax."
  (interactive)
  (beginning-of-line)
  (if (bobp)
      (indent-line-to 0)		   ; First line indented to column 0
    (let ((not-indented t) cur-indent)
      (if (or (looking-at cp-decr-indent-fold-t-re)
	      (let ((case-fold-search nil)) (looking-at cp-decr-indent-re))) ; If the line we are looking at is the end of a block, then decrease the indentation
	  (progn
	    (save-excursion
	      (forward-line -1)
	      (setq cur-indent (- (current-indentation) qe-indent)))
	    (if (< cur-indent 0) ; We can't indent past the left margin
		(setq cur-indent 0)))
	(save-excursion
	  (while  not-indented ; Iterate backwards until we find an indentation hint
	    (forward-line -1)
	    (if (looking-at cp-deindent-fold-t-re) ; This hint indicates that we need to indent at the level of the "/" token
		(progn
		  (setq cur-indent (current-indentation))
		  (setq not-indented nil))
	      (if (or (looking-at cp-indent-fold-t-re)
		      (let ((case-fold-search nil)) (looking-at cp-indent-re))) ; This hint indicates that we need to indent an extra level
		  (progn
		    (setq cur-indent (+ (current-indentation) qe-indent)) ; Do the actual indenting
		    (setq not-indented nil))
		(if (bobp)
		    (setq not-indented nil)))))))
      (if cur-indent
	  (indent-line-to cur-indent)
	(indent-line-to 0))))) ; If we didn't see an indentation hint, then allow no indentation


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the cp-mode as derived-mode

(define-derived-mode cp-mode prog-mode 
  "QE-cp.x"  
  "Major mode for editing Qunatum-ESPRESSO input files (QE-cp.x mode)"  
  (setq font-lock-defaults '((cp-keywords) nil t))
  (set (make-local-variable 'indent-line-function) 'cp-indent-line)
  
  ;; define the syntax of comments
  (setq comment-start "!")
  (setq comment-end "")
  (modify-syntax-entry ?!  "< b" cp-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" cp-mode-syntax-table)
  (modify-syntax-entry ?=  " " cp-mode-syntax-table) ;; treat "=" non symbol constituent
  ;; end
  )

;; free memory

(setq cp-namelists nil)
(setq cp-vars nil)
(setq cp-cards nil)
(setq cp-flags nil)
(setq cp-open-supercards   nil)
(setq cp-closed-supercards nil)


(require 'qe-funcs)

(provide 'cp-mode)

