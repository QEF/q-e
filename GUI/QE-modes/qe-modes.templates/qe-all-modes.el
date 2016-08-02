;;
;; generic QE mode: knows all QE namlists, cards, etc..
;;
(autoload 'qe-mode' "qe-mode.el"
  "Major mode for editing Quantum ESPRESSO input files" t)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; specific modes: each knows only its own namlists, cards, etc..
;;

(autoload 'pw-mode' "pw-mode.el"
  "Major mode for editing Quantum ESPRESSO pw.x input files" t)

(autoload 'ph-mode' "ph-mode.el"
  "Major mode for editing Quantum ESPRESSO ph.x input files" t)

(autoload 'pp-mode' "pp-mode.el"
  "Major mode for editing Quantum ESPRESSO pp.x input files" t)

(autoload 'projwfc-mode' "projwfc-mode.el"
  "Major mode for editing Quantum ESPRESSO projwfc.x input files" t)

(autoload 'atomic-mode' "atomic-mode.el"
  "Major mode for editing Quantum ESPRESSO ld1.x input files" t)

(autoload 'cp-mode' "cp-mode.el"
  "Major mode for editing Quantum ESPRESSO cp.x input files" t)


(provide 'qe-all-modes)
