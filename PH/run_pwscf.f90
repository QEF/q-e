!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE run_pwscf(do_band)
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the pwscf program called from the
  ! ... phonon code.
  !
  !
  USE control_flags,   ONLY : conv_ions, twfcollect
  USE basis,           ONLY : starting_wfc, starting_pot, startingconfig
  USE io_files,        ONLY : prefix, tmp_dir
  USE lsda_mod,        ONLY : nspin
  USE input_parameters,ONLY : pseudo_dir
  USE control_flags,   ONLY : restart
  USE qpoint,          ONLY : xq
  USE control_ph,      ONLY : done_bands, reduce_io, recover, tmp_dir_phq, &
                              ext_restart, bands_computed, newgrid
  USE save_ph,         ONLY : tmp_dir_save
  !
  USE acfdtest,      ONLY : acfdt_is_active, acfdt_num_der, ir_point, delta_vrs
  USE scf,           ONLY : vrs

 !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: do_band
  !
  CHARACTER(LEN=256) :: dirname, file_base_in, file_base_out
  LOGICAL :: exst
  !
  CALL start_clock( 'PWSCF' )
  !
  CALL clean_pw( .FALSE. )
  !
  CALL close_files(.true.)
  !
  ! From now on, work only on the _ph virtual directory
  !
  tmp_dir=tmp_dir_phq
  !
  ! ... Setting the values for the nscf run
  !
  startingconfig    = 'input'
  starting_pot      = 'file'
  starting_wfc      = 'atomic'
  restart = ext_restart
  pseudo_dir= TRIM( tmp_dir_save ) // TRIM( prefix ) // '.save'
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
     done_bands=.TRUE.
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
  CALL stop_clock( 'PWSCF' )
  !
  RETURN
END SUBROUTINE run_pwscf
