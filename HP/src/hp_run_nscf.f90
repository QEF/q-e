!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE hp_run_nscf (do_band)
  !-----------------------------------------------------------------------
  !
  ! This is the main driver of the PWscf program called from the HP code.
  !
  USE control_flags,   ONLY : conv_ions, restart, iverbosity
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : prefix, tmp_dir, wfc_dir, seqopn
  USE lsda_mod,        ONLY : nspin
  USE check_stop,      ONLY : check_stop_now
  USE fft_types,       ONLY : fft_type_allocate
  USE fft_base,        ONLY : dffts, dfftp
  USE cell_base,       ONLY : at, bg
  USE gvect,           ONLY : gcutm
  USE gvecs,           ONLY : gcutms
  USE io_global,       ONLY : stdout
  USE scf,             ONLY : vrs
  USE mp_bands,        ONLY : intra_bgrp_comm, nyfft
  USE qpoint,          ONLY : xq
  USE control_lr,      ONLY : lgamma
  USE lr_symm_base,    ONLY : nsymq, invsymq
  USE ldaU_hp,         ONLY : tmp_dir_save, tmp_dir_hp
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: do_band
  LOGICAL :: exst
  INTEGER :: verbosity_save
  !
  CALL start_clock( 'hp_run_nscf' )
  !
  CALL clean_pw (.FALSE.)
  !
  CALL close_files (.TRUE.)
  !
  ! From now on, work only on the HP directory
  !
  wfc_dir = tmp_dir_hp
  tmp_dir = tmp_dir_hp
  !
  ! Setting the values for the NSCF run
  !
  startingconfig = 'input'
  starting_pot   = 'file'
  starting_wfc   = 'atomic'
  restart        = .FALSE.
  conv_ions      = .TRUE.
  !
  ! iverbosity is used by the PWscf routines
  IF (iverbosity.LE.2) THEN
     ! temporarily change the value of iverbosity
     ! in order to have less output from the PWscf routines
     verbosity_save = iverbosity
     iverbosity = 0
  ENDIF
  !
  IF (lgamma) THEN
     WRITE(stdout,'(/5x,"Performing NSCF calculation at all points k...")')
  ELSE
     WRITE(stdout,'(/5x,"Performing NSCF calculation at all points k and k+q...")')
  ENDIF
  !
  CALL fft_type_allocate ( dfftp, at, bg, gcutm,  intra_bgrp_comm, nyfft=nyfft )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft)
  !
  CALL setup_nscf ( .FALSE., xq, .FALSE. )
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
  WRITE(stdout,'(5x,"Done!")')
  !
  CALL stop_clock( 'hp_run_nscf' )
  !
  RETURN
  ! 
END SUBROUTINE hp_run_nscf
