!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE kcw_run_nscf (do_band)
  !-----------------------------------------------------------------------
  !
  ! This is the main driver of the PWscf program called from the KCW code.
  !
  USE control_flags,   ONLY : conv_ions, restart, iverbosity, isolve
  USE starting_scf,    ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : tmp_dir, wfc_dir
  USE fft_types,       ONLY : fft_type_allocate
  USE fft_base,        ONLY : dffts, dfftp
  USE cell_base,       ONLY : at, bg
  USE gvect,           ONLY : gcutm
  USE gvecs,           ONLY : gcutms
  USE mp_bands,        ONLY : intra_bgrp_comm, nyfft
  USE qpoint,          ONLY : xq
  USE control_lr,      ONLY : ethr_nscf
  USE control_kcw,     ONLY : tmp_dir_kcwq, kcw_iverbosity
  USE klist,           ONLY : nelec
  USE control_kcw,     ONLY : irr_bz
  USE control_kcw,     ONLY : mp1, mp2, mp3
  USE klist,           ONLY : xk, wk
  USE start_k,         ONLY : init_start_k, nks_start, xk_start
  USE start_k,         ONLY : nk1, nk2, nk3, k1, k2, k3
  USE io_global,       ONLY : stdout
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: do_band
  INTEGER :: verbosity_save, ik
  !
  CALL start_clock( 'kcw_run_nscf' )
  !
  CALL clean_pw (.FALSE.)
  !
  CALL close_files (.TRUE.)
  !
  ! From now on, work only on the KC directory
  !
  wfc_dir = tmp_dir_kcwq
  tmp_dir = tmp_dir_kcwq
  !
  ! Setting the values for the NSCF run
  !
  startingconfig = 'input'
  starting_pot   = 'file'
  starting_wfc   = 'atomic'
  restart        = .FALSE.
  conv_ions      = .TRUE.
  ethr_nscf      = 1.0D-9 / nelec 
  isolve         = 0
  !
  ! iverbosity is used by the PWscf routines
  IF (iverbosity.LE.2) THEN
     ! temporarily change the value of iverbosity
     ! in order to have less output from the PWscf routines
     verbosity_save = iverbosity
     iverbosity = 0
  ENDIF
  !
  CALL fft_type_allocate ( dfftp, at, bg, gcutm,  intra_bgrp_comm, nyfft=nyfft )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft)
  !
  ! we need the IBZ of k if sym.
  !
  IF( irr_bz) THEN
    !WRITE(*,*) "nks_start: ", nks_start
    !
    ! need to setup nk1, nk2, nk3 from start_k
    !if we don't do it, the code will crash when reading the charge density.
    !
    !nk1 = mp1
    !nk2 = mp2
    !nk3 = mp3
    !k1 = 0.D0
    !k2 = 0.D0
    !k3 = 0.D0
    !
   CALL init_start_k( mp1, mp2, mp3, 0, 0, 0, 'cartesian', &
   mp1*mp2*mp3, xk, wk ) 
   
   IF (kcw_iverbosity .gt. 2) THEN 
     WRITE(stdout ,'(8X, "SYM :nks_start =", i5)') nks_start
     !DO ik=1, nks_start 
      WRITE(stdout,'(8X, "     xk_start(",i5, " )=", 3F12.4)') (ik, xk_start(:,ik), ik=1, nks_start)
     !END DO
   ENDIF
   CALL setup_nscf ( .TRUE., xq, .FALSE. ) 
   ! the second false restrict to IBZ(q) the k points
  ELSE 
    CALL setup_nscf ( .FALSE., xq, .TRUE. ) 
  END IF
  !
  CALL init_run()
  !
  IF (do_band) THEN
     CALL non_scf_ph()
     CALL punch( 'all' )
  ENDIF
  !
  IF (iverbosity.EQ.0) iverbosity = verbosity_save 
  !
  CALL close_files(.TRUE.)
  !
  CALL stop_clock( 'kcw_run_nscf' )
  !
  RETURN
  ! 
END SUBROUTINE kcw_run_nscf
