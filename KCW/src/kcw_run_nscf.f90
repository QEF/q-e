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
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : tmp_dir, wfc_dir
  USE fft_types,       ONLY : fft_type_allocate
  USE fft_base,        ONLY : dffts, dfftp
  USE cell_base,       ONLY : at, bg
  USE gvect,           ONLY : gcutm
  USE gvecs,           ONLY : gcutms
  USE mp_bands,        ONLY : intra_bgrp_comm, nyfft
  USE qpoint,          ONLY : xq
  USE control_lr,      ONLY : ethr_nscf
  USE control_kcw,     ONLY : tmp_dir_kcwq
  USE klist,           ONLY : nelec
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: do_band
  INTEGER :: verbosity_save
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
  CALL setup_nscf ( .FALSE., xq, .TRUE. )
  !
  CALL init_run()
  !
  IF (do_band) THEN
     CALL non_scf()
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
