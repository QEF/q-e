!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE recover_mod

  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE

  INTEGER ::  iunrec=99

  PUBLIC :: write_rec, read_rec, clean_recover

CONTAINS

  !-----------------------------------------------------------------------
  SUBROUTINE write_rec(where, irr, dr2, iter, convt, npe, dvscfin, &
       drhoscfh, dbecsum)
    !-----------------------------------------------------------------------
    !
    !  This routine saves the information needed to recover the phonon
    !
    USE kinds, ONLY : DP
    USE ions_base, ONLY : nat
    USE uspp_param, ONLY : nhm
    USE lsda_mod,  ONLY : nspin
    USE units_ph, ONLY : this_pcxpsi_is_on_file
    USE noncollin_module, ONLY : nspin_mag
    USE fft_base, ONLY : dfftp
    USE uspp, ONLY : okvan, nlcc_any
    USE phus, ONLY : int1, int2
    USE control_ph, ONLY : where_rec, rec_code, reduce_io, current_iq
    USE ph_restart, ONLY : ph_writefile
    USE efield_mod, ONLY : zstareu0, zstarue0
    USE io_files, ONLY : seqopn

    USE lrus,   ONLY : int3
    USE eqv,    ONLY : drhoscfs
    USE qpoint, ONLY : nksq

    IMPLICIT NONE
    CHARACTER(LEN=10), INTENT(IN) :: where
    INTEGER, INTENT(IN) :: irr, iter, npe
    LOGICAL, INTENT(IN) :: convt
    REAL(DP), INTENT(IN) :: dr2
    COMPLEX(DP), INTENT(IN) :: dvscfin(dfftp%nnr,nspin_mag,npe)
    COMPLEX(DP), INTENT(IN), OPTIONAL :: drhoscfh (dfftp%nnr, nspin_mag, npe)
    COMPLEX(DP), INTENT(IN), OPTIONAL :: dbecsum((nhm*(nhm+1))/2,nat,nspin_mag,npe)

    INTEGER :: ierr
    LOGICAL :: exst
    CALL start_clock ('write_rec')
    where_rec=where
    CALL ph_writefile('status_ph',current_iq,0,ierr)
    IF (where=='done_drhod') CALL ph_writefile('data_dyn',current_iq,irr,ierr)
    CALL seqopn (iunrec, 'recover', 'unformatted', exst)
    !
    ! info on current iteration (iter=0 potential mixing not available)
    !
    IF (reduce_io.or.convt) THEN
       WRITE (iunrec) 0, dr2, convt
    ELSE
       WRITE (iunrec) iter, dr2, convt
    ENDIF
    WRITE (iunrec) this_pcxpsi_is_on_file
    WRITE (iunrec) zstareu0, zstarue0
    WRITE (iunrec) dvscfin
    IF (PRESENT(drhoscfh).AND.convt.AND.nlcc_any) WRITE (iunrec) drhoscfh
    IF (convt.AND.ALLOCATED(drhoscfs)) WRITE(iunrec) drhoscfs
    IF (PRESENT(dbecsum)) WRITE(iunrec) dbecsum
    IF (okvan) WRITE (iunrec) int1, int2, int3

    CLOSE (UNIT = iunrec, STATUS = 'keep')

    rec_code = 0
    CALL stop_clock ('write_rec')

    RETURN
  END SUBROUTINE write_rec

  SUBROUTINE read_rec(dr2, iter0, npe, dvscfin, dvscfins, drhoscfh, dbecsum)
    !
    !  General restart reading routine
    !
    USE kinds, ONLY : DP
    USE ions_base, ONLY : nat
    USE uspp_param, ONLY : nhm
    USE gvecs, ONLY : doublegrid
    USE fft_base, ONLY : dfftp, dffts
    USE uspp,  ONLY : okvan, nlcc_any
    USE lsda_mod, ONLY : nspin
    USE noncollin_module, ONLY : noncolin, nspin_mag
    USE units_ph, ONLY : this_pcxpsi_is_on_file
    USE control_ph, ONLY : ext_recover, convt
    USE efield_mod, ONLY : zstareu0, zstarue0
    USE phus, ONLY : int1, int2
    USE io_files, ONLY : seqopn

    USE lrus, ONLY : int3
    USE eqv,  ONLY : drhoscfs

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: iter0
    INTEGER, INTENT(IN)  :: npe
    REAL(DP), INTENT(OUT) :: dr2
    COMPLEX(DP), INTENT(OUT) :: dvscfin (dfftp%nnr, nspin_mag, npe)
    COMPLEX(DP), INTENT(OUT) :: dvscfins (dffts%nnr, nspin_mag, npe)
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: drhoscfh (dfftp%nnr, nspin_mag, npe)
    COMPLEX(DP), INTENT(OUT), OPTIONAL :: dbecsum((nhm*(nhm+1))/2,nat,nspin_mag,npe)

    INTEGER :: is, ipol
    LOGICAL :: exst

    CALL start_clock ('read_rec')
    CALL seqopn (iunrec, 'recover', 'unformatted', exst)
    READ (iunrec) iter0, dr2, convt
    READ (iunrec) this_pcxpsi_is_on_file
    READ (iunrec) zstareu0, zstarue0
    READ (iunrec) dvscfin
    IF (convt.AND.nlcc_any) READ(iunrec) drhoscfh
    IF (convt.AND.ALLOCATED(drhoscfs)) READ(iunrec) drhoscfs
    IF (PRESENT(dbecsum)) READ(iunrec) dbecsum
    IF (okvan) THEN
       READ (iunrec) int1, int2, int3
       IF (noncolin) THEN
          CALL set_int12_nc(0)
          CALL set_int3_nc(npe)
       END IF
    END IF
    CLOSE (UNIT = iunrec, STATUS = 'keep')
    IF (doublegrid) THEN
       DO is=1,nspin_mag
          DO ipol=1,npe
             CALL cinterpolate (dvscfin(1,is,ipol), dvscfins(1,is,ipol), -1)
          END DO
       END DO
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
