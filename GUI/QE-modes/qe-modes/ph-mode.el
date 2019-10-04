;; ph-mode.el
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

;; This is the `ph-mode', a major mode for composing the Quantum ESPRESSO
;; QE-ph.x input files. For the installation and usage, see the
;; user_guide.pdf in the Doc/ subdirectory of the original package
;; (quick installation instructions are also available in the README
;; file of the original package).

;;; Code:

(require 'font-lock)
(require 'regexp-opt)

(defvar ph-mode-hook nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; basic variables

;; ph's supercards (if any)
(defvar ph-open-supercards   (list ))
(defvar ph-closed-supercards (list ))
  
;; ph's namelists
(defvar ph-namelists (list "&INPUTPH" ))
(defvar qe-end-namelist (list "&END" "/"))

;; ph's variables
(defvar ph-vars (list "alpha_mix" "amass" "asr" "dek" "drho_star" "dvscf_star" "electron_phonon" "elop" "epsil" "eth_ns" "eth_rps" "fildrho" "fildvscf" "fildyn" "fpol" "k1" "k2" "k3" "last_irr" "last_q" "ldiag" "ldisp" "lnoloc" "low_directory_check" "lqdir" "lraman" "lrpa" "lshift_q" "max_seconds" "modenum" "nat_todo" "niter_ph" "nk1" "nk2" "nk3" "nmix_ph" "nogg" "nq1" "nq2" "nq3" "only_init" "outdir" "prefix" "q2d" "q_in_band_form" "qplot" "read_dns_bare" "recover" "reduce_io" "search_sym" "start_irr" "start_q" "tr2_ph" "trans" "verbosity" "zeu" "zue" ))

;; ph's cards & keywords
(defvar ph-cards (list "__NO-CARDS" ))

;; ph's flags
(defvar ph-flags (list ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; derived variables
  
(defvar ph-open-supercards-regexp   (regexp-opt ph-open-supercards   'symbols)) ; may not exists
(defvar ph-closed-supercards-regexp (regexp-opt ph-closed-supercards 'symbols)) ; may not exists

(defvar ph-cards-regexp (regexp-opt
			    (append ph-cards ph-open-supercards) 'symbols))
(defvar ph-flags-regexp (regexp-opt ph-flags 'symbols))

(defvar ph-namelist-face (cons (regexp-opt (append ph-namelists qe-end-namelist) 'symbols) font-lock-function-name-face))
(defvar ph-variable-face (cons (regexp-opt ph-vars 'symbols) font-lock-variable-name-face))

;; logical values as constants
(defvar qe-logic-face (cons (regexp-opt (list ".t." ".true." ".f." ".false.")) font-lock-constant-face))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; regexp for indentation
(defvar ph-decr-indent-fold-t-re (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
(defvar ph-decr-indent-re        (concat "^[ \t]*" (regexp-opt
						       (append ph-cards ph-open-supercards ph-closed-supercards) t)))
;;
(defvar ph-deindent-fold-t-re    (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
;;
(defvar ph-indent-fold-t-re      (concat "^[ \t]*" (regexp-opt ph-namelists t)))
(defvar ph-indent-re             (concat "^[ \t]*" (regexp-opt ph-cards     t)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; supercards, cards and flags are case sensitive -- here are the corresponding matchers

(defun ph-closed-supercards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward ph-closed-supercards-regexp limit 'no-error)))

(defun ph-cards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward ph-cards-regexp limit 'no-error)))

(defun ph-flags-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward ph-flags-regexp limit 'no-error)))


(font-lock-add-keywords 'ph-mode (list
				     ph-namelist-face 
				     ph-variable-face
				     qe-logic-face
				     '("," . font-lock-builtin-face)
				     '("(" . font-lock-builtin-face)
				     '(")" . font-lock-builtin-face)
				     ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; register the keywords

(font-lock-add-keywords 'ph-mode '(
				      (ph-closed-supercards-matcher 1 font-lock-preprocessor-face t)
				      (ph-cards-matcher 1 font-lock-keyword-face t)
				      (ph-flags-matcher 1 font-lock-type-face    t)
				      ))

;;(defvar ph-keywords '(ph-namelist-face ph-variable-face))
(defvar ph-keywords '(((list "") . font-lock-constant-face)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the sytnax of strings

(defvar ph-mode-syntax-table
  (let ((table (make-syntax-table)))
    (modify-syntax-entry ?\' "\"'"  table)
    (modify-syntax-entry ?\" "\"\"" table)
    table)
  "Syntax table in use in `ph-mode' buffers.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; code for auto-indenting

(defvar qe-indent 3)
(defun ph-indent-line ()
  "Indent current line according to ph input syntax."
  (interactive)
  (beginning-of-line)
  (if (bobp)
      (indent-line-to 0)		   ; First line indented to column 0
    (let ((not-indented t) cur-indent)
      (if (or (looking-at ph-decr-indent-fold-t-re)
	      (let ((case-fold-search nil)) (looking-at ph-decr-indent-re))) ; If the line we are looking at is the end of a block, then decrease the indentation
	  (progn
	    (save-excursion
	      (forward-line -1)
	      (setq cur-indent (- (current-indentation) qe-indent)))
	    (if (< cur-indent 0) ; We can't indent past the left margin
		(setq cur-indent 0)))
	(save-excursion
	  (while  not-indented ; Iterate backwards until we find an indentation hint
	    (forward-line -1)
	    (if (looking-at ph-deindent-fold-t-re) ; This hint indicates that we need to indent at the level of the "/" token
		(progn
		  (setq cur-indent (current-indentation))
		  (setq not-indented nil))
	      (if (or (looking-at ph-indent-fold-t-re)
		      (let ((case-fold-search nil)) (looking-at ph-indent-re))) ; This hint indicates that we need to indent an extra level
		  (progn
		    (setq cur-indent (+ (current-indentation) qe-indent)) ; Do the actual indenting
		    (setq not-indented nil))
		(if (bobp)
		    (setq not-indented nil)))))))
      (if cur-indent
	  (indent-line-to cur-indent)
	(indent-line-to 0))))) ; If we didn't see an indentation hint, then allow no indentation


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the ph-mode as derived-mode

(define-derived-mode ph-mode prog-mode 
  "QE-ph.x"  
  "Major mode for editing Qunatum-ESPRESSO input files (QE-ph.x mode)"  
  (setq font-lock-defaults '((ph-keywords) nil t))
  (set (make-local-variable 'indent-line-function) 'ph-indent-line)
  
  ;; define the syntax of comments
  (setq comment-start "!")
  (setq comment-end "")
  (modify-syntax-entry ?!  "< b" ph-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" ph-mode-syntax-table)
  (modify-syntax-entry ?=  " " ph-mode-syntax-table) ;; treat "=" non symbol constituent
  ;; end
  )

;; free memory

(setq ph-namelists nil)
(setq ph-vars nil)
(setq ph-cards nil)
(setq ph-flags nil)
(setq ph-open-supercards   nil)
(setq ph-closed-supercards nil)


(require 'qe-funcs)

(provide 'ph-mode)

