!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------
SUBROUTINE lr_init_nfo()
  !-------------------------------------------------------------
  !
  !  This subroutine prepares several variables which are needed in the
  !  TDDFPT program:
  !  1) Initialization of ikks, ikqs, and nksq.
  !  2) EELS: Calculate phases associated with a q vector.
  !  3) Compute the number of occupied bands for each k point.
  !  4) Computes alpha_pv (needed by orthogonalize.f90 when lgauss=.true.)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, tau
  USE klist,                ONLY : nks,xk,ngk,igk_k
  USE wvfct,                ONLY : nbnd
  USE realus,               ONLY : real_space
  USE lr_variables,         ONLY : lr_verbosity, eels, restart, &
                                   size_evc, tmp_dir_lr
  USE io_global,            ONLY : stdout
  USE constants,            ONLY : tpi, eps8
  USE noncollin_module,     ONLY : npol
  USE gvect,                ONLY : ngm, g
  USE cell_base,            ONLY : at, bg, omega, tpiba2
  USE lsda_mod,             ONLY : current_spin, nspin
  USE wvfct,                ONLY : npwx, wg
  USE gvecw,                ONLY : gcutw
  USE io_files,             ONLY : prefix, iunwfc, nwordwfc, wfc_dir
  USE gvecs,                ONLY : doublegrid
  USE fft_base,             ONLY : dfftp 
  USE uspp,                 ONLY : vkb, okvan, nkb
  USE wavefunctions_module, ONLY : evc
  USE becmod,               ONLY : calbec, allocate_bec_type
  USE lrus,                 ONLY : becp1
  USE control_lr,           ONLY : alpha_pv
  USE qpoint,               ONLY : xq, ikks, ikqs, nksq, eigqts
  USE eqv,                  ONLY : evq
  USE buffers,              ONLY : open_buffer, get_buffer
  USE control_flags,        ONLY : io_level
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  REAL(kind=DP) :: arg
  INTEGER :: na,  & ! dummy index which runs over atoms 
             ik,  & ! dummy index which runs over k points
             ikk, & ! index of the point k
             npw    ! number of the plane-waves at point k
  LOGICAL :: exst, exst_mem 
  ! logical variable to check file exists
  ! logical variable to check file exists in memory
  !
  ! 1) Initialization of ikks, ikqs, and nksq.
  !    Inspired by PH/initialize_ph.f90.
  !
  IF (eels) THEN
     !
     ! EELS (q!=0)
     !
     ! nksq is the number of k-points, NOT including k+q points
     !
     nksq = nks / 2
     !
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik = 1, nksq
        ikks(ik) = 2 * ik - 1
        ikqs(ik) = 2 * ik
     ENDDO
     !
  ELSE
     !
     ! Optical case (q=0)
     !
     nksq = nks
     !
     ! Note: in this case the arrays ikks and ikqs 
     ! are needed only in h_prec and ch_psi_all.
     !
     ALLOCATE(ikks(nksq), ikqs(nksq))
     DO ik = 1, nksq
        ikks(ik) = ik
        ikqs(ik) = ik
     ENDDO
     !
  ENDIF
  !
  ! 2) EELS-specific operations
  !
  IF (eels) THEN
     !
     size_evc = nbnd * npwx * npol * nksq
     nwordwfc = nbnd * npwx * npol
     !
     ! Open file to read the wavefunctions at k and k+q points 
     ! after the nscf calculation.
     !
     IF (restart) wfc_dir = tmp_dir_lr
     !
     CALL open_buffer (iunwfc, 'wfc', nwordwfc, io_level, exst_mem, exst)
     ! 
     IF (.NOT.exst .AND. .NOT.exst_mem) THEN
        CALL errore ('lr_init_nfo', 'file '//trim(prefix)//'.wfc not found', 1)
     ENDIF
     !
     ! If restart=.true. recalculate the small group of q.
     !
     IF (restart) CALL lr_smallgq (xq)
     !
     ! USPP-specific initializations
     !
     IF (okvan) THEN
        !
        ! Calculate phases associated with the q vector
        !
        ALLOCATE (eigqts(nat))
        !
        DO na = 1, nat
           arg = ( xq(1) * tau(1,na) + &
                   xq(2) * tau(2,na) + &
                   xq(3) * tau(3,na) ) * tpi
           eigqts(na) = cmplx(cos(arg), -sin(arg) ,kind=DP)
        ENDDO
        !
        ! Calculate becp1 = <vkb|evc>
        !
        ALLOCATE (becp1(nksq))
        !
        DO ik = 1, nksq
           !
           CALL allocate_bec_type (nkb,nbnd,becp1(ik))
           !
           ikk = ikks(ik)
           npw = ngk(ikk)
           !
           ! Read the wavefunction evc
           CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
           !
           ! Calculate beta-functions vkb at point k
           CALL init_us_2(npw, igk_k(1,ikk), xk(1,ikk), vkb)
           !
           ! Calculate becp1=<vkb|evc>
           CALL calbec (npw, vkb, evc, becp1(ik))
           !
        ENDDO
        !
     ENDIF
     !
  ENDIF
  !
  ! 3) Compute the number of occupied bands for each k point
  !
  CALL setup_nbnd_occ()
  !
  ! 4) Compute alpha_pv
  !
  IF (eels) THEN
     !
     alpha_pv = 0.0d0
     !
  ELSE
     !
     CALL setup_alpha_pv()
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE lr_init_nfo
