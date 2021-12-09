!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE lr_run_nscf( )
  !--------------------------------------------------------------------------
  !
  ! This is the main driver of the PWscf program called from the TDDFPT code.
  ! Inspired by PH/run_nscf.f90
  !
  ! Created by Iurii Timrov (2013)
  !
  USE control_flags,   ONLY : conv_ions, restart
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : prefix, tmp_dir, wfc_dir, seqopn
  USE io_global,       ONLY : stdout
  USE lr_variables,    ONLY : tmp_dir_lr
  USE mp_bands,        ONLY : ntask_groups, intra_bgrp_comm, nyfft
  USE fft_types,       ONLY : fft_type_allocate
  USE fft_base,        ONLY : dffts, dfftp
  USE cell_base,       ONLY : at, bg
  USE gvect,           ONLY : gcutm
  USE gvecs,           ONLY : gcutms
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst
  !
  CALL start_clock( 'lr_run_nscf' )
  !
  WRITE( stdout, '(/5x,28("-")," NSCF calculation ",28("-"))')
  !   
  CALL clean_pw(.FALSE.)
  !   
  CALL close_files(.TRUE.)
  !
  ! From now on, work only on the _ph directory
  !
  wfc_dir = tmp_dir_lr
  tmp_dir = tmp_dir_lr
  !
  ! Setting the values for the nscf run
  !
  startingconfig    = 'input'
  starting_pot      = 'file'
  starting_wfc      = 'atomic'
  !
  ! Do not confuse the below restart (from control_flags) with
  ! the one from the TDDFPT code (from lr_variables).
  !
  restart = .false.
  conv_ions = .true.  ! IT: maybe this is not needed
  !
  CALL fft_type_allocate ( dfftp, at, bg, gcutm,  intra_bgrp_comm, nyfft=nyfft )
  CALL fft_type_allocate ( dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft )
  !
  ! Initialize variables for the non-scf calculations at k
  ! and k+q required by the linear response calculation at finite q.
  !
  CALL lr_setup_nscf ()
  !
  ! The variables igk_k and ngk are re-set up here (because there 
  ! are not only poins k but also points k+q) through the path:
  ! init_run -> allocate_wfc_k -> init_igk
  ! 
  CALL init_run()
  !
  ! Non-scf calculation
  !
  CALL non_scf()
  !
  ! Save information for further processing
  !
  CALL punch( 'all' ) 
  !
  CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !
  CALL close_files(.TRUE.)
  !   
  CALL stop_clock( 'lr_run_nscf' )
  !
  RETURN
  !
END SUBROUTINE lr_run_nscf
