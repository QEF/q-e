!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
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
  USE control_flags,   ONLY : conv_ions, twfcollect
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : prefix, tmp_dir, wfc_dir, seqopn
  USE lsda_mod,        ONLY : nspin
  USE control_flags,   ONLY : restart
  USE control_ph,      ONLY : reduce_io, tmp_dir_phq, ext_restart
  USE save_ph,         ONLY : tmp_dir_save
  USE io_global,       ONLY : stdout
  USE fft_base,        ONLY : dffts
  USE mp_bands,        ONLY : ntask_groups
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst
  !
  CALL start_clock( 'lr_run_nscf' )
  !   
  CALL clean_pw(.FALSE.)
  !   
  CALL close_files(.TRUE.)
  !
  ! From now on, work only on the _ph directory
  !
  wfc_dir = tmp_dir_phq
  tmp_dir = tmp_dir_phq
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
  ! Initialize variables for the non-scf calculations at k
  ! and k+q required by the linear response calculation at finite q.
  !
  CALL lr_setup_nscf ()
  ! 
  CALL init_run()
  !
  ! Non-scf calculation
  !
  CALL non_scf()
  !
  IF (.not.reduce_io) THEN
     !
     twfcollect = .FALSE.
     CALL punch( 'all' ) 
     !
  ENDIF
  !
  CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  !
  CALL close_files(.TRUE.)
  !   
  !  PWscf has run with task groups if available, but in the TDDFPT (phonon) 
  !  they are used only in some places. In that case it is activated.
  !
  IF (ntask_groups > 1) dffts%have_task_groups = .FALSE.
  !
  CALL stop_clock( 'lr_run_nscf' )
  !
  RETURN
  !
END SUBROUTINE lr_run_nscf
