!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE recover_mod
  !
  !! Module for phonon recovery.
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE

  INTEGER ::  iunrec=99

  PUBLIC :: write_rec, read_rec, clean_recover

CONTAINS

  !-----------------------------------------------------------------------
  SUBROUTINE write_rec( where, irr, dr2, iter, convt, dfpt_data )
    !-----------------------------------------------------------------------
    !! This routine saves the information needed to recover the phonon.
    !
    USE kinds, ONLY : DP
    USE units_ph, ONLY : this_pcxpsi_is_on_file
    USE uspp, ONLY : okvan
    USE paw_variables, ONLY : okpaw
    USE phus, ONLY : int1, int2
    USE control_lr, ONLY : where_rec, rec_code, reduce_io
    USE control_ph, ONLY : current_iq
    USE ph_restart, ONLY : ph_writefile
    USE efield_mod, ONLY : zstareu0, zstarue0
    USE io_files, ONLY : seqopn
    USE dfpt_type, ONLY : dfpt_data_type
    USE lrus,   ONLY : int3

    IMPLICIT NONE
    CHARACTER(LEN=10), INTENT(IN) :: where
    INTEGER, INTENT(IN) :: irr, iter
    LOGICAL, INTENT(IN) :: convt
    REAL(DP), INTENT(IN) :: dr2
    TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data

    INTEGER :: ierr
    LOGICAL :: exst

    !
    ! Write dynmat.X.Y if computed, regardless io_level.
    ! Dynamical matrices should be written also when reduce_io is true, to recover 
    ! from accomplished runs on images. Otherwise, dynmat.X.Y files are not found 
    ! and trans calculation recovers from scratch.
    !
    CALL start_clock ('write_rec')
    where_rec=where
    CALL ph_writefile('status_ph',current_iq,0,ierr)
    IF (where=='done_drhod') CALL ph_writefile('data_dyn',current_iq,irr,ierr)
    IF (reduce_io) THEN
       CALL stop_clock ('write_rec')   
       RETURN
    ENDIF
    CALL seqopn (iunrec, 'recover', 'unformatted', exst)
    !
    ! info on current iteration (iter=0 potential mixing not available)
    !
    IF (convt) THEN
       WRITE (iunrec) 0, dr2, convt
    ELSE
       WRITE (iunrec) iter, dr2, convt
    ENDIF
    WRITE (iunrec) this_pcxpsi_is_on_file
    WRITE (iunrec) zstareu0, zstarue0
    WRITE (iunrec) dfpt_data%dvscfp
    WRITE (iunrec) dfpt_data%drhop
    WRITE (iunrec) dfpt_data%drhos
    IF (okpaw) WRITE (iunrec) dfpt_data%dbecsum
    IF (okvan) WRITE (iunrec) int1, int2, int3

    CLOSE (UNIT = iunrec, STATUS = 'keep')

    CALL stop_clock ('write_rec')

    RETURN
  END SUBROUTINE write_rec
  
  !-----------------------------------------------------------------------------------
  SUBROUTINE read_rec( dr2, iter0, dfpt_data)
    !--------------------------------------------------------------------------------
    !! General restart reading routine.
    !
    USE kinds, ONLY : DP
    USE gvecs, ONLY : doublegrid
    USE fft_base, ONLY : dfftp, dffts
    USE fft_interfaces, ONLY : fft_interpolate
    USE uspp,  ONLY : okvan
    USE paw_variables, ONLY : okpaw
    USE noncollin_module, ONLY : noncolin, nspin_mag
    USE units_ph, ONLY : this_pcxpsi_is_on_file
    USE control_lr, ONLY : convt
    USE control_ph, ONLY : ext_recover
    USE efield_mod, ONLY : zstareu0, zstarue0
    USE phus, ONLY : int1, int2
    USE io_files, ONLY : seqopn
    USE dfpt_type, ONLY : dfpt_data_type
    USE lrus, ONLY : int3

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: iter0
    REAL(DP), INTENT(OUT) :: dr2
    TYPE(dfpt_data_type), INTENT(INOUT) :: dfpt_data

    INTEGER :: is, ipol
    LOGICAL :: exst

    CALL start_clock ('read_rec')
    CALL seqopn (iunrec, 'recover', 'unformatted', exst)
    READ (iunrec) iter0, dr2, convt
    READ (iunrec) this_pcxpsi_is_on_file
    READ (iunrec) zstareu0, zstarue0
    READ (iunrec) dfpt_data%dvscfp
    READ (iunrec) dfpt_data%drhop
    READ (iunrec) dfpt_data%drhos
    !
    ! If not PAW, dbecsum is used only for the calculation of drhop, so we do not need it
    ! for restart. For PAW, we need dbecsum for restart.
    IF (okpaw) READ (iunrec) dfpt_data%dbecsum
    !
    IF (okvan) THEN
       READ (iunrec) int1, int2, int3
       IF (noncolin) THEN
          CALL set_int12_nc(0)
          CALL set_int3_nc(dfpt_data%npert)
       END IF
    END IF
    CLOSE (UNIT = iunrec, STATUS = 'keep')
    !
    IF (doublegrid) THEN
       DO is=1,nspin_mag
          DO ipol=1,dfpt_data%npert
             CALL fft_interpolate (dfftp, dfpt_data%dvscfp(:,is,ipol), dffts, dfpt_data%dvscfs(:,is,ipol))
          END DO
       END DO
    ELSE
       CALL zcopy(dfftp%nnr*nspin_mag*dfpt_data%npert, dfpt_data%dvscfp, 1, dfpt_data%dvscfs, 1)
    END IF
    ext_recover=.FALSE.
    CALL stop_clock ('read_rec')

    RETURN
  END SUBROUTINE read_rec

  SUBROUTINE clean_recover()
    !
    USE io_files, ONLY : seqopn
    !
    IMPLICIT NONE
    LOGICAL :: exst
    !
    CALL seqopn( iunrec, 'recover', 'UNFORMATTED', exst )
    !
    CLOSE( UNIT = iunrec, STATUS = 'DELETE' )
    !
  END SUBROUTINE clean_recover

END MODULE recover_mod
