$header

;;; Commentary:

;; This is the `$mode-mode', a major mode for composing the $modeName scripts.
;; For the installation and usage, see the user_guide.pdf in the Doc/
;; subdirectory of the original package (quick installation
;; instructions are also available in the README file of the original
;; package).

;;; Code:

(require 'font-lock)
(require 'regexp-opt)

(defvar ${mode}-mode-hook nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; basic variables

;; QE's supercards (if any)
(defvar qe-open-supercards   (list $open_supercards))
(defvar qe-closed-supercards (list $closed_supercards))
  
;; QE's namelists
(defvar qe-namelists (list $namelists))

;; QE's variables
(defvar qe-vars (list $vars))

;; QE's cards & keywords
(defvar qe-cards (list $cards))

;; QE's flags
(defvar qe-flags (list $flags))

;; ${modeName}'s cmds
(defvar $mode-cmds (list $cmds))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; derived variables
  
(defvar qe-open-supercards-regexp   (regexp-opt qe-open-supercards   'symbols)) ; may not exists
(defvar qe-closed-supercards-regexp (regexp-opt qe-closed-supercards 'symbols)) ; may not exists

(defvar qe-cards-regexp (regexp-opt
			    (append qe-cards qe-open-supercards) 'symbols))
(defvar qe-flags-regexp (regexp-opt qe-flags 'symbols))

(defvar qe-namelist-face (cons (regexp-opt (append qe-namelists) 'symbols) font-lock-function-name-face))
(defvar $mode-cmds-face (cons (regexp-opt (append $mode-cmds) 'symbols) font-lock-function-name-face))
(defvar qe-variable-face (cons (regexp-opt qe-vars 'symbols) font-lock-variable-name-face))

;; logical values as constants
(defvar qe-logic-face (cons (regexp-opt (list ".t." ".true." ".f." ".false.")) font-lock-constant-face))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; supercards, cards and flags are case sensitive -- here are the corresponding matchers

(defun qe-closed-supercards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward qe-closed-supercards-regexp limit 'no-error)))

(defun qe-cards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward qe-cards-regexp limit 'no-error)))

(defun qe-flags-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward qe-flags-regexp limit 'no-error)))


(font-lock-add-keywords '$mode-mode (list
				     qe-namelist-face 
				     $mode-cmds-face 
				     qe-variable-face
				     qe-logic-face
				     '("," . font-lock-builtin-face)
				     '("(" . font-lock-builtin-face)
				     '(")" . font-lock-builtin-face)
				     ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; register the keywords

(font-lock-add-keywords '$mode-mode '(
				      (qe-closed-supercards-matcher 1 font-lock-preprocessor-face t)
                                      (qe-cards-matcher 1 font-lock-function-name-face t)
				      (qe-flags-matcher 1 font-lock-type-face    t)
				      ))

;;(defvar qe-keywords '(qe-namelist-face qe-variable-face))
(defvar qe-keywords '(((list "") . font-lock-constant-face)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the sytnax of strings

(defvar $mode-mode-syntax-table
  (let ((table (make-syntax-table)))
    (modify-syntax-entry ?\' "\"'"  table)
    (modify-syntax-entry ?\" "\"\"" table)
    table)
  "Syntax table in use in `$mode-mode' buffers.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the $mode-mode as derived-mode

(define-derived-mode $mode-mode tcl-mode 
  "$modeName-mode"  
  "Major mode for editing $modeName scripts"    
  )

;; free memory

(setq qe-namelists nil)
(setq pwtk-cmds nil)
(setq qe-vars nil)
(setq qe-cards nil)
(setq qe-flags nil)
(setq qe-open-supercards   nil)
(setq qe-closed-supercards nil)


(require 'qe-funcs)
(provide '$mode-mode)
