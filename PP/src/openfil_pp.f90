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
  USE basis,          ONLY : natomwfc
  USE ldaU,           ONLY : nwfcU
  USE control_flags,  ONLY : twfcollect
  USE io_files,       ONLY : prefix, iunwfc, diropn, &
                             nwordwfc, nwordatwfc, nwordwfcU
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
  ! NOTE: in post-processing codes, wavefunctions are still opened using
  !       "diropn" and not "open_buffer" (there is no real advantage in
  !       using buffers; there would be one if wavefunctions in collected
  !       format were read into a buffer instead of being written to file
  !       in distributed format, but this is not yet done). In order to have
  !       a uniform definition, nwordwfc is defined as in pwscf as the number
  !       of COMPLEX WORDS of the wavefunction packet.
  !
  nwordwfc = nbnd * npwx * npol
  nwordwfcU = nwfcU * npwx * npol
  nwordatwfc = natomwfc * npwx * npol
  !
  CALL diropn( iunwfc, 'wfc', 2*nwordwfc, exst )
  IF ( .not. exst ) &
     CALL errore ('openfil_pp','file '//trim( prefix )//'.wfc'//' not found',1)
  !
  RETURN
  !
END SUBROUTINE openfil_pp
