;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Purpose:
;;          this is the snippet for the $HOME/.emacs file. Edit it
;;          according to your needs and insert the respective content
;;          into $HOME/.emacs file
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;
;; QE modes: BEGIN



(add-to-list 'load-path (expand-file-name "~/.emacs.d/qe-modes/"))
(require 'qe-modes)


;; automatic filename recognition
;; (BEWARE: more general patterns must be specified first)

;; automatically open the *.in files with generic QE mode
(add-to-list 'auto-mode-alist '("\\.in\\'" . qe-mode))

;;  automatically open the the pw*.in, scf*.in, nscf*in, relax*in,
;;  vc-relax*.in, md*.in, vc-md*.in files by pw.x mode
(add-to-list 'auto-mode-alist
	     '("\\(pw\\|n?scf\\|\\(?:vc-\\)?\\(?:md\\|relax\\)\\).*\\.in\\'" . pw-mode))

;; automatically open the neb*.in files with neb.x mode
(add-to-list 'auto-mode-alist '("neb.*\\.in\\'" . neb-mode))

;; automatically open the cp*.in files with cp.x mode
(add-to-list 'auto-mode-alist '("cp.*\\.in\\'" . ph-mode))

;; automatically open the ph*.in files with ph.x mode
(add-to-list 'auto-mode-alist '("ph.*\\.in\\'" . ph-mode))

;; automatically open the ld1*.in files with ld1 mode
(add-to-list 'auto-mode-alist '("ld1.*\\.in\\'" . ld1-mode))

;; automatically open the pp*.in files with pp.x mode
(add-to-list 'auto-mode-alist '("pp.*\\.in\\'" . pp-mode))


;;; default indentation offset is 3; uncomment below line and set the
;;; value accordingly if you prefer other value
;
;(setq qe-indent 4)


;;; uncomment below lines to disable the auto-indentation ...
;;; ( are you really sure you want to do this ? )
;
;;(dolist (hook '(qe-mode-hook
;;                pw-mode-hook
;;                cp-mode-hook
;;                ph-mode-hook
;;                ld1-mode-hook
;;		  pp-mode-hook))
;;  (add-hook hook (lambda () (setq indent-line-function 'indent-relative))))


;; QE modes: END
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
