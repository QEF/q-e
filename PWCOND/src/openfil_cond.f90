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
  USE klist,            ONLY : nks
  USE io_files,         ONLY : prefix, iunpun, iunat, iunsat, iunwfc, iunigk, &
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
  !
  nwordwfc = nbnd*npwx*npol
  !
  ! ... iunwfc= 10: read/write wfc from/to file
  ! ... iunwfc=-10: copy wfc to/from RAM
  !
  IF ( io_level > 0 ) THEN
     iunwfc = 10
  ELSE
     iunwfc =-10
  END IF
  CALL open_buffer( iunwfc, 'wfc', nwordwfc, nks, exst )
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
  USE control_flags,    ONLY : io_level
  !
  IMPLICIT NONE
  !
  !
  IF ( io_level > 0 ) THEN
     iunwfc = 10
  ELSE
     iunwfc =-10
  END IF
  CALL close_buffer( iunwfc, 'keep' )
  !
  RETURN
  !
END SUBROUTINE closefil_cond



