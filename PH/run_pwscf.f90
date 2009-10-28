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
  USE control_ph,      ONLY : done_bands, reduce_io, recover, tmp_dir_ph
  USE save_ph,         ONLY : tmp_dir_save
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: dirname, file_base_in, file_base_out
  !
  LOGICAL, INTENT(IN) :: do_band
  !
  LOGICAL :: exst
  !
  CALL start_clock( 'PWSCF' )
  !
  CALL clean_pw( .FALSE. )
  !
  CALL close_files()
  !
  ! From now on, work only on the _ph virtual directory
  !
  tmp_dir=tmp_dir_ph
  !
  ! ... Setting the values for the nscf run
  !
  startingconfig    = 'input'
  starting_pot      = 'file'
  starting_wfc      = 'atomic'
  restart = recover
  pseudo_dir= TRIM( tmp_dir_save ) // TRIM( prefix ) // '.save'
  CALL restart_from_file()
  conv_ions=.true.
  !
  CALL setup_nscf (xq)
  CALL init_run()
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
  !
  CALL close_files()
  !
  CALL stop_clock( 'PWSCF' )
  !
  RETURN
END SUBROUTINE run_pwscf
