$header

;;; Commentary:

;; This is the `$mode-mode', a major mode for composing the Quantum ESPRESSO
;; $modeName input files. For the installation and usage, see the
;; user_guide.pdf in the Doc/ subdirectory of the original package
;; (quick installation instructions are also available in the README
;; file of the original package).

;;; Code:

(require 'font-lock)
(require 'regexp-opt)

(defvar ${mode}-mode-hook nil)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; basic variables

;; ${mode}'s supercards (if any)
(defvar $mode-open-supercards   (list $open_supercards))
(defvar $mode-closed-supercards (list $closed_supercards))
  
;; ${mode}'s namelists
(defvar $mode-namelists (list $namelists))
(defvar qe-end-namelist (list "&END" "/"))

;; ${mode}'s variables
(defvar $mode-vars (list $vars))

;; ${mode}'s cards & keywords
(defvar $mode-cards (list $cards))

;; ${mode}'s flags
(defvar $mode-flags (list $flags))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; derived variables
  
(defvar $mode-open-supercards-regexp   (regexp-opt $mode-open-supercards   'symbols)) ; may not exists
(defvar $mode-closed-supercards-regexp (regexp-opt $mode-closed-supercards 'symbols)) ; may not exists

(defvar $mode-cards-regexp (regexp-opt
			    (append $mode-cards $mode-open-supercards) 'symbols))
(defvar $mode-flags-regexp (regexp-opt $mode-flags 'symbols))

(defvar $mode-namelist-face (cons (regexp-opt (append $mode-namelists qe-end-namelist) 'symbols) font-lock-function-name-face))
(defvar $mode-variable-face (cons (regexp-opt $mode-vars 'symbols) font-lock-variable-name-face))

;; logical values as constants
(defvar qe-logic-face (cons (regexp-opt (list ".t." ".true." ".f." ".false.")) font-lock-constant-face))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; regexp for indentation
(defvar $mode-decr-indent-fold-t-re (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
(defvar $mode-decr-indent-re        (concat "^[ \t]*" (regexp-opt
						       (append $mode-cards $mode-open-supercards $mode-closed-supercards) t)))
;;
(defvar $mode-deindent-fold-t-re    (concat "^[ \t]*" (regexp-opt qe-end-namelist t)))
;;
(defvar $mode-indent-fold-t-re      (concat "^[ \t]*" (regexp-opt $mode-namelists t)))
(defvar $mode-indent-re             (concat "^[ \t]*" (regexp-opt $mode-cards     t)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; supercards, cards and flags are case sensitive -- here are the corresponding matchers

(defun $mode-closed-supercards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward $mode-closed-supercards-regexp limit 'no-error)))

(defun $mode-cards-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward $mode-cards-regexp limit 'no-error)))

(defun $mode-flags-matcher (limit)
  (let ((case-fold-search nil))
    (re-search-forward $mode-flags-regexp limit 'no-error)))


(font-lock-add-keywords '$mode-mode (list
				     $mode-namelist-face 
				     $mode-variable-face
				     qe-logic-face
				     '("," . font-lock-builtin-face)
				     '("(" . font-lock-builtin-face)
				     '(")" . font-lock-builtin-face)
				     ))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; register the keywords

(font-lock-add-keywords '$mode-mode '(
				      ($mode-closed-supercards-matcher 1 font-lock-preprocessor-face t)
				      ($mode-cards-matcher 1 font-lock-keyword-face t)
				      ($mode-flags-matcher 1 font-lock-type-face    t)
				      ))

;;(defvar $mode-keywords '($mode-namelist-face $mode-variable-face))
(defvar $mode-keywords '(((list "") . font-lock-constant-face)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the sytnax of strings

(defvar $mode-mode-syntax-table
  (let ((table (make-syntax-table)))
    (modify-syntax-entry ?\' "\"'"  table)
    (modify-syntax-entry ?\" "\"\"" table)
    table)
  "Syntax table in use in `$mode-mode' buffers.")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; code for auto-indenting

(defvar qe-indent 3)
(defun $mode-indent-line ()
  "Indent current line according to $mode input syntax."
  (interactive)
  (beginning-of-line)
  (if (bobp)
      (indent-line-to 0)		   ; First line indented to column 0
    (let ((not-indented t) cur-indent)
      (if (or (looking-at $mode-decr-indent-fold-t-re)
	      (let ((case-fold-search nil)) (looking-at $mode-decr-indent-re))) ; If the line we are looking at is the end of a block, then decrease the indentation
	  (progn
	    (save-excursion
	      (forward-line -1)
	      (setq cur-indent (- (current-indentation) qe-indent)))
	    (if (< cur-indent 0) ; We can't indent past the left margin
		(setq cur-indent 0)))
	(save-excursion
	  (while  not-indented ; Iterate backwards until we find an indentation hint
	    (forward-line -1)
	    (if (looking-at $mode-deindent-fold-t-re) ; This hint indicates that we need to indent at the level of the "/" token
		(progn
		  (setq cur-indent (current-indentation))
		  (setq not-indented nil))
	      (if (or (looking-at $mode-indent-fold-t-re)
		      (let ((case-fold-search nil)) (looking-at $mode-indent-re))) ; This hint indicates that we need to indent an extra level
		  (progn
		    (setq cur-indent (+ (current-indentation) qe-indent)) ; Do the actual indenting
		    (setq not-indented nil))
		(if (bobp)
		    (setq not-indented nil)))))))
      (if cur-indent
	  (indent-line-to cur-indent)
	(indent-line-to 0))))) ; If we didn't see an indentation hint, then allow no indentation


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; define the $mode-mode as derived-mode

(define-derived-mode $mode-mode prog-mode 
  "$modeName"  
  "Major mode for editing Qunatum-ESPRESSO input files ($modeName mode)"  
  (setq font-lock-defaults '(($mode-keywords) nil t))
  (set (make-local-variable 'indent-line-function) '$mode-indent-line)
  
  ;; define the syntax of comments
  (setq comment-start "!")
  (setq comment-end "")
  (modify-syntax-entry ?!  "< b" $mode-mode-syntax-table)
  (modify-syntax-entry ?\n "> b" $mode-mode-syntax-table)
  (modify-syntax-entry ?=  " " $mode-mode-syntax-table) ;; treat "=" non symbol constituent
  ;; end
  )

;; free memory

(setq $mode-namelists nil)
(setq $mode-vars nil)
(setq $mode-cards nil)
(setq $mode-flags nil)
(setq $mode-open-supercards   nil)
(setq $mode-closed-supercards nil)


(require 'qe-funcs)

(provide '$mode-mode)
