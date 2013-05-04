!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE openfil_pp()
  !----------------------------------------------------------------------------
  !
  ! ... This routine opens all files needed to the self consistent run,
  ! ... sets various file names, units, record lengths
  !
  USE kinds,          ONLY : DP
  USE wvfct,          ONLY : nbnd, npwx
  USE control_flags,  ONLY:  twfcollect
  USE io_files,       ONLY : prefix, iunwfc, nwordwfc, diropn
  USE noncollin_module, ONLY : npol
  !
  IMPLICIT NONE
  !
  LOGICAL       :: exst
  !
  twfcollect=.false.
  !
  ! ... nwordwfc is the record length for the direct-access file
  ! ... containing wavefunctions
  !
  ! FIXME: in post-processing codes, wavefunctions are still opened using
  !        "diropn" and not "open_buffer" (there is no real advantage in
  !        using buffers; there would be one if wavefunctions in collected
  !        format were read into a buffer instead of being written to file
  !        in distributed format, but this is not yet done) As a consequence
  !        nwordwfc is twice as big as in pwscf. This is not a problem as
  !        long as PP does not attempt to use use "get_buffer" instead of
  !        "davcio" (e.g. via routines in PW). PG May 2013
  !
  nwordwfc = 2 * nbnd * npwx * npol
  CALL diropn( iunwfc, 'wfc', nwordwfc, exst )
  !
  IF ( .not. exst ) THEN
     CALL errore ('openfil_pp','file '//trim( prefix )//'.wfc'//' not found',1)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE openfil_pp
