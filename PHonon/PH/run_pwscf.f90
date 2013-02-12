!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE run_pwscf(do_band, iq)
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the pwscf program called from the
  ! ... phonon code.
  !
  !
  USE control_flags,   ONLY : conv_ions, twfcollect
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : prefix, tmp_dir, wfc_dir, seqopn
  USE lsda_mod,        ONLY : nspin
  USE control_flags,   ONLY : restart
  USE fft_base,        ONLY : dffts
  USE modes,           ONLY : minus_q, nsymq, invsymq
  USE disp,            ONLY : lgamma_iq
  USE qpoint,          ONLY : xq
  USE control_ph,      ONLY : reduce_io, recover, tmp_dir_phq, &
                              ext_restart, bands_computed, newgrid
  USE save_ph,         ONLY : tmp_dir_save
  !
  USE grid_irr_iq,     ONLY : done_bands
  USE acfdtest,        ONLY : acfdt_is_active, acfdt_num_der, ir_point, delta_vrs
  USE scf,             ONLY : vrs
  USE mp_global,       ONLY : get_ntask_groups

 !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: do_band
  INTEGER, INTENT(IN) :: iq
  !
  LOGICAL :: exst
  !
  CALL start_clock( 'PWSCF' )
  !
  IF (done_bands(iq)) THEN
     CALL clean_pw( .TRUE. )
     CALL close_files(.true.)
     wfc_dir=tmp_dir_phq
     tmp_dir=tmp_dir_phq
     CALL read_file()
!
!   This routine is here only for invsymq. It seems that this flag is not
!   used. If it is useful to somebody please tell me or the call to this routine
!   will be removed.
!
     IF (.NOT.lgamma_iq(iq)) CALL set_small_group_of_q(nsymq,invsymq,minus_q)
     RETURN
  ENDIF
  !
  CALL clean_pw( .FALSE. )
  !
  CALL close_files(.true.)
  !
  ! From now on, work only on the _ph virtual directory
  !
  wfc_dir=tmp_dir_phq
  tmp_dir=tmp_dir_phq
  ! ... Setting the values for the nscf run
  !
  startingconfig    = 'input'
  starting_pot      = 'file'
  starting_wfc      = 'atomic'
  restart = ext_restart
  CALL restart_from_file()
  conv_ions=.true.
  !
  CALL setup_nscf ( newgrid, xq )
  CALL init_run()
!!!!!!!!!!!!!!!!!!!!!!!! ACFDT TEST !!!!!!!!!!!!!!!!
  IF (acfdt_is_active) THEN
    ! ACFDT mumerical derivative test: modify the potential
    IF (acfdt_num_der) vrs(ir_point,1)=vrs(ir_point,1) + delta_vrs
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!END OF ACFDT TEST !!!!!!!!!!!!!!!!
!
  IF (do_band) CALL electrons()
  !
  IF (.NOT.reduce_io.and.do_band) THEN
     twfcollect=.FALSE.
     CALL punch( 'all' )
  ENDIF
  !
  CALL seqopn( 4, 'restart', 'UNFORMATTED', exst )
  CLOSE( UNIT = 4, STATUS = 'DELETE' )
  ext_restart=.FALSE.
  !
  CALL close_files(.true.)
  !

  bands_computed=.TRUE.
!
!  PWscf has run with task groups if available, but in the phonon 
!  they are not used, apart in particular points. In that case it is
!  activated.
!
  IF (get_ntask_groups()>1) dffts%have_task_groups=.FALSE.
  !
  CALL stop_clock( 'PWSCF' )
  !
  RETURN
END SUBROUTINE run_pwscf
