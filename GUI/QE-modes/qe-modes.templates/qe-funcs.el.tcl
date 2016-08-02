$header

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


$utility_functions


$keyword_functions


(provide 'qe-funcs)


