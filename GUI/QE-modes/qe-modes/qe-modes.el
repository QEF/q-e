;; qe-modes.el
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

(autoload 'pw-mode' "pw-mode.elc"
  "Major mode for editing Quantum ESPRESSO pw.x input files" t)

(autoload 'cp-mode' "cp-mode.elc"
  "Major mode for editing Quantum ESPRESSO cp.x input files" t)

(autoload 'pp-mode' "pp-mode.elc"
  "Major mode for editing Quantum ESPRESSO pp.x input files" t)

(autoload 'ld1-mode' "ld1-mode.elc"
  "Major mode for editing Quantum ESPRESSO ld1.x input files" t)

(autoload 'neb-mode' "neb-mode.elc"
  "Major mode for editing Quantum ESPRESSO neb.x input files" t)

(autoload 'ph-mode' "ph-mode.elc"
  "Major mode for editing Quantum ESPRESSO ph.x input files" t)



(provide 'qe-modes)

