;; pp-mode.el
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

;; This is the `pp-mode', a major mode for composing the Quantum ESPRESSO
;; QE-pp.x input files. For the installation and usage, see the
;; user_guide.pdf in the Doc/ subdirectory of the original package
;; (quick installation instructions are also available in the README
;; file of the original package).

;;; Code:

(require 'font-lock)
(require 'regexp-opt)

(defvar pp-mode-hook nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; basic variables

;; pp's supercards (if any)
(defvar pp-open-supercards   (list ))
(defvar pp-closed-supercards (list ))
  
;; pp's namelists
(defvar pp-namelists (list "&INPUTPP" "&PLOT" ))
(defvar qe-end-namelist (list "&END" "/"))

;; pp's variables
(defvar pp-vars (list "degauss_ldos" "delta_e" "e1" "e2" "e3" "emax" "emin" "fileout" "filepp" "filplot" "iflag" "interpolation" "kband" "kpoint" "lsign" "nfile" "nx" "ny" "nz" "outdir" "output_format" "plot_num" "prefix" "radius" "sample_bias" "spin_component" "weight" "x0" ))

;; pp's cards & keywords
(defvar pp-cards (list "__NO-CARDS" ))

;; pp's flags
(defvar pp-flags (list ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; derived variables
  
(defvar pp-open-supercards-regexp   (regexp-opt pp-open-supercards   'symbols)) ; may not exists
(defvar pp-closed-supercards-regexp (regexp-opt pp-closed-supercards 'symbols)) ; may not exists

(defvar pp-cards-regexp (regexp-opt
			    (append pp-cards pp-open-supercards) 'symbols))
(defvar pp-flags-regexp (regexp-opt pp-flags 'symbols))

(defvar pp-namelist-face (cons (regexp-opt (append pp-namelists qe-end-namelist) 'symbols) font-lock-function-name-face))
(defvar pp-variable-face (cons (regexp-opt pp-vars 'symbols) font-lock-variable-name-face))

;; logical values as constants
(defvar qe-logic-face (cons (regexp-opt (list ".t." ".true." ".f." ".false.")) font-lock-constant-face))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; regexp for indentation
(defvar pp-decr-indent-fold-t-re (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
(defvar pp-decr-indent-re        (concat "^[ \t]*" (regexp-opt
						       (append pp-cards pp-open-supercards pp-closed-supercards) t)))
;;
(defvar pp-deindent-fold-t-re    (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
;;
(defvar pp-indent-fold-t-re      (concat "^[ \t]*" (regexp-opt pp-namelists t)))
(defvar pp-indent-re             (concat "^[ \t]*" (regexp-opt pp-cards     t)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; supercards, cards and flags are case sensitive -- here are the corresponding matchers

(defun pp-closed-supercards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward pp-closed-supercards-regexp limit 'no-error)))

(defun pp-cards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward pp-cards-regexp limit 'no-error)))

(defun pp-flags-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward pp-flags-regexp limit 'no-error)))


(font-lock-add-keywords 'pp-mode (list
				     pp-namelist-face 
				     pp-variable-face
				     qe-logic-face
				     '("," . font-lock-builtin-face)
				     '("(" . font-lock-builtin-face)
				     '(")" . font-lock-builtin-face)
				     ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; register the keywords

(font-lock-add-keywords 'pp-mode '(
				      (pp-closed-supercards-matcher 1 font-lock-preprocessor-face t)
				      (pp-cards-matcher 1 font-lock-keyword-face t)
				      (pp-flags-matcher 1 font-lock-type-face    t)
				      ))

;;(defvar pp-keywords '(pp-namelist-face pp-variable-face))
(defvar pp-keywords '(((list "") . font-lock-constant-face)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the sytnax of strings

(defvar pp-mode-syntax-table
  (let ((table (make-syntax-table)))
    (modify-syntax-entry ?\' "\"'"  table)
    (modify-syntax-entry ?\" "\"\"" table)
    table)
  "Syntax table in use in `pp-mode' buffers.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; code for auto-indenting

(defvar qe-indent 3)
(defun pp-indent-line ()
  "Indent current line according to pp input syntax."
  (interactive)
  (beginning-of-line)
  (if (bobp)
      (indent-line-to 0)		   ; First line indented to column 0
    (let ((not-indented t) cur-indent)
      (if (or (looking-at pp-decr-indent-fold-t-re)
	      (let ((case-fold-search nil)) (looking-at pp-decr-indent-re))) ; If the line we are looking at is the end of a block, then decrease the indentation
	  (progn
	    (save-excursion
	      (forward-line -1)
	      (setq cur-indent (- (current-indentation) qe-indent)))
	    (if (< cur-indent 0) ; We can't indent past the left margin
		(setq cur-indent 0)))
	(save-excursion
	  (while  not-indented ; Iterate backwards until we find an indentation hint
	    (forward-line -1)
	    (if (looking-at pp-deindent-fold-t-re) ; This hint indicates that we need to indent at the level of the "/" token
		(progn
		  (setq cur-indent (current-indentation))
		  (setq not-indented nil))
	      (if (or (looking-at pp-indent-fold-t-re)
		      (let ((case-fold-search nil)) (looking-at pp-indent-re))) ; This hint indicates that we need to indent an extra level
		  (progn
		    (setq cur-indent (+ (current-indentation) qe-indent)) ; Do the actual indenting
		    (setq not-indented nil))
		(if (bobp)
		    (setq not-indented nil)))))))
      (if cur-indent
	  (indent-line-to cur-indent)
	(indent-line-to 0))))) ; If we didn't see an indentation hint, then allow no indentation


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the pp-mode as derived-mode

(define-derived-mode pp-mode prog-mode 
  "QE-pp.x"  
  "Major mode for editing Qunatum-ESPRESSO input files (QE-pp.x mode)"  
  (setq font-lock-defaults '((pp-keywords) nil t))
  (set (make-local-variable 'indent-line-function) 'pp-indent-line)
  
  ;; define the syntax of comments
  (setq comment-start "!")
  (setq comment-end "")
  (modify-syntax-entry ?!  "< b" pp-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" pp-mode-syntax-table)
  (modify-syntax-entry ?=  " " pp-mode-syntax-table) ;; treat "=" non symbol constituent
  ;; end
  )

;; free memory

(setq pp-namelists nil)
(setq pp-vars nil)
(setq pp-cards nil)
(setq pp-flags nil)
(setq pp-open-supercards   nil)
(setq pp-closed-supercards nil)


(require 'qe-funcs)

(provide 'pp-mode)

