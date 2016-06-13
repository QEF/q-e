!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE openfil_cond()
  !----------------------------------------------------------------------------
  !
  ! ... This routine opens some files needed by pwcond,
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE wvfct,            ONLY : nbnd, npwx
  USE io_files,         ONLY : prefix, iunpun, iunsat, iunwfc, &
                               nwordwfc, nwordatwfc, iunefield, &
                               iunefieldm, iunefieldp
  USE noncollin_module, ONLY : npol
  USE mp_global,        ONLY : kunit
  USE buffers,          ONLY : open_buffer
  USE control_flags,    ONLY : io_level
  !
  IMPLICIT NONE
  !
  LOGICAL            :: exst
  !
  ! ... nwordwfc is the record length (IN COMPLEX WORDS)
  ! ... for the direct-access file containing wavefunctions
  ! ... io_level > 0 : open a file; io_level <= 0 : open a buffer
  !
  nwordwfc = nbnd*npwx*npol
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, io_level, exst )
  !
  RETURN
  !
END SUBROUTINE openfil_cond

!----------------------------------------------------------------------------
SUBROUTINE closefil_cond()
  !----------------------------------------------------------------------------
  !
  ! ... This routine close the files opened by pwcond
  !
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : iunwfc
  USE buffers,          ONLY : close_buffer
  !
  IMPLICIT NONE
  !
  CALL close_buffer( iunwfc, 'keep' )
  !
  RETURN
  !
END SUBROUTINE closefil_cond



