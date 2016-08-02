$header

;;; Commentary:

;; This is the master `qe-modes.el' file, which loads various Quantum
;; ESPRESSO major modes. Add the following into your user-init-file
;; (e.g. ~/.emacs):
;;
;;    (add-to-list 'load-path "/directory/of/your/qe-modes/")
;;    (require 'qe-modes)

;;; Code:

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; load generic Quantum ESPRESSO mode: it knows all QE namelists, cards, etc.
;;
(autoload 'qe-mode' "qe-mode.elc"
  "Major mode for editing Quantum ESPRESSO input files" t)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; load specific modes: each knows only its own namelists, cards, etc.
;;

$autoload_specific_modes

(provide 'qe-modes)
