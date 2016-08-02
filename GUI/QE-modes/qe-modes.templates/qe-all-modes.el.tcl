$header


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; load generic Quantum ESPRESSO mode: it knows all QE namelists, cards, etc.
;;
(autoload 'qe-mode' "qe-mode.el"
  "Major mode for editing Quantum ESPRESSO input files" t)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; load specific modes: each knows only its own namelists, cards, etc.
;;

$autoload_specific_modes

(provide 'qe-modes)
