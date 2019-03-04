;; ld1-mode.el
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

;; This is the `ld1-mode', a major mode for composing the Quantum ESPRESSO
;; QE-ld1.x (atomic) input files. For the installation and usage, see the
;; user_guide.pdf in the Doc/ subdirectory of the original package
;; (quick installation instructions are also available in the README
;; file of the original package).

;;; Code:

(require 'font-lock)
(require 'regexp-opt)

(defvar ld1-mode-hook nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; basic variables

;; ld1's supercards (if any)
(defvar ld1-open-supercards   (list ))
(defvar ld1-closed-supercards (list ))
  
;; ld1's namelists
(defvar ld1-namelists (list "&INPUT" "&INPUTP" "&TEST" ))
(defvar qe-end-namelist (list "&END" "/"))

;; ld1's variables
(defvar ld1-vars (list "atom" "author" "beta" "cau_fact" "config" "configts" "decut" "deld" "dft" "dx" "ecutmax" "ecutmin" "emaxld" "eminld" "file_beta" "file_charge" "file_chi" "file_core" "file_pseudo" "file_pseudopw" "file_qvan" "file_recon" "file_screen" "file_wfcaegen" "file_wfcncgen" "file_wfcusgen" "frozen_core" "isic" "iswitch" "latt" "lgipaw_reconstruction" "lloc" "lpaw" "lsave_wfc" "lsd" "lsdts" "lsmall" "max_out_wfc" "nconf" "new_core_ps" "nlcc" "nld" "noscf" "prefix" "pseudotype" "rcloc" "rcore" "rcutv" "rel" "rel_dist" "relpert" "rho0" "rlderiv" "rm" "rmatch_augfun" "rmatch_augfun_nc" "rmax" "rpwe" "rytoev_fact" "title" "tm" "tr2" "use_paw_as_gipaw" "vdw" "verbosity" "which_augfun" "write_coulomb" "xmin" "zed" "zval" ))

;; ld1's cards & keywords
(defvar ld1-cards (list "__NO-CARDS" ))

;; ld1's flags
(defvar ld1-flags (list ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; derived variables
  
(defvar ld1-open-supercards-regexp   (regexp-opt ld1-open-supercards   'symbols)) ; may not exists
(defvar ld1-closed-supercards-regexp (regexp-opt ld1-closed-supercards 'symbols)) ; may not exists

(defvar ld1-cards-regexp (regexp-opt
			    (append ld1-cards ld1-open-supercards) 'symbols))
(defvar ld1-flags-regexp (regexp-opt ld1-flags 'symbols))

(defvar ld1-namelist-face (cons (regexp-opt (append ld1-namelists qe-end-namelist) 'symbols) font-lock-function-name-face))
(defvar ld1-variable-face (cons (regexp-opt ld1-vars 'symbols) font-lock-variable-name-face))

;; logical values as constants
(defvar qe-logic-face (cons (regexp-opt (list ".t." ".true." ".f." ".false.")) font-lock-constant-face))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; regexp for indentation
(defvar ld1-decr-indent-fold-t-re (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
(defvar ld1-decr-indent-re        (concat "^[ \t]*" (regexp-opt
						       (append ld1-cards ld1-open-supercards ld1-closed-supercards) t)))
;;
(defvar ld1-deindent-fold-t-re    (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
;;
(defvar ld1-indent-fold-t-re      (concat "^[ \t]*" (regexp-opt ld1-namelists t)))
(defvar ld1-indent-re             (concat "^[ \t]*" (regexp-opt ld1-cards     t)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; supercards, cards and flags are case sensitive -- here are the corresponding matchers

(defun ld1-closed-supercards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward ld1-closed-supercards-regexp limit 'no-error)))

(defun ld1-cards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward ld1-cards-regexp limit 'no-error)))

(defun ld1-flags-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward ld1-flags-regexp limit 'no-error)))


(font-lock-add-keywords 'ld1-mode (list
				     ld1-namelist-face 
				     ld1-variable-face
				     qe-logic-face
				     '("," . font-lock-builtin-face)
				     '("(" . font-lock-builtin-face)
				     '(")" . font-lock-builtin-face)
				     ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; register the keywords

(font-lock-add-keywords 'ld1-mode '(
				      (ld1-closed-supercards-matcher 1 font-lock-preprocessor-face t)
				      (ld1-cards-matcher 1 font-lock-keyword-face t)
				      (ld1-flags-matcher 1 font-lock-type-face    t)
				      ))

;;(defvar ld1-keywords '(ld1-namelist-face ld1-variable-face))
(defvar ld1-keywords '(((list "") . font-lock-constant-face)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the sytnax of strings

(defvar ld1-mode-syntax-table
  (let ((table (make-syntax-table)))
    (modify-syntax-entry ?\' "\"'"  table)
    (modify-syntax-entry ?\" "\"\"" table)
    table)
  "Syntax table in use in `ld1-mode' buffers.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; code for auto-indenting

(defvar qe-indent 3)
(defun ld1-indent-line ()
  "Indent current line according to ld1 input syntax."
  (interactive)
  (beginning-of-line)
  (if (bobp)
      (indent-line-to 0)		   ; First line indented to column 0
    (let ((not-indented t) cur-indent)
      (if (or (looking-at ld1-decr-indent-fold-t-re)
	      (let ((case-fold-search nil)) (looking-at ld1-decr-indent-re))) ; If the line we are looking at is the end of a block, then decrease the indentation
	  (progn
	    (save-excursion
	      (forward-line -1)
	      (setq cur-indent (- (current-indentation) qe-indent)))
	    (if (< cur-indent 0) ; We can't indent past the left margin
		(setq cur-indent 0)))
	(save-excursion
	  (while  not-indented ; Iterate backwards until we find an indentation hint
	    (forward-line -1)
	    (if (looking-at ld1-deindent-fold-t-re) ; This hint indicates that we need to indent at the level of the "/" token
		(progn
		  (setq cur-indent (current-indentation))
		  (setq not-indented nil))
	      (if (or (looking-at ld1-indent-fold-t-re)
		      (let ((case-fold-search nil)) (looking-at ld1-indent-re))) ; This hint indicates that we need to indent an extra level
		  (progn
		    (setq cur-indent (+ (current-indentation) qe-indent)) ; Do the actual indenting
		    (setq not-indented nil))
		(if (bobp)
		    (setq not-indented nil)))))))
      (if cur-indent
	  (indent-line-to cur-indent)
	(indent-line-to 0))))) ; If we didn't see an indentation hint, then allow no indentation


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the ld1-mode as derived-mode

(define-derived-mode ld1-mode prog-mode 
  "QE-ld1.x (atomic)"  
  "Major mode for editing Qunatum-ESPRESSO input files (QE-ld1.x (atomic) mode)"  
  (setq font-lock-defaults '((ld1-keywords) nil t))
  (set (make-local-variable 'indent-line-function) 'ld1-indent-line)
  
  ;; define the syntax of comments
  (setq comment-start "!")
  (setq comment-end "")
  (modify-syntax-entry ?!  "< b" ld1-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" ld1-mode-syntax-table)
  (modify-syntax-entry ?=  " " ld1-mode-syntax-table) ;; treat "=" non symbol constituent
  ;; end
  )

;; free memory

(setq ld1-namelists nil)
(setq ld1-vars nil)
(setq ld1-cards nil)
(setq ld1-flags nil)
(setq ld1-open-supercards   nil)
(setq ld1-closed-supercards nil)


(require 'qe-funcs)

(provide 'ld1-mode)

