!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE close_files()
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes for a new scf calculation.
  !
  USE control_flags, ONLY : order
  USE io_files,      ONLY : prefix, iunwfc, iunoldwfc, iunoldwfc2, iunigk
  USE mp_global,     ONLY : intra_image_comm
  USE mp,            ONLY : mp_barrier
  !
  IMPLICIT NONE
  !
  LOGICAL :: flag
  LOGICAL :: exst
  !
  !  ... iunwfc contains wavefunctions and is kept open during
  !  ... the execution - close and save the file
  !
  CLOSE( UNIT = iunwfc, STATUS = 'KEEP' )
  !
  ! ... iunigk is kept open during the execution - close and remove
  !
  CLOSE( UNIT = iunigk, STATUS = 'DELETE' )
  !
  CALL mp_barrier( intra_image_comm )  
  !
#ifdef __T3E
  !
  ! ... set streambuffers off
  !
  CALL set_d_stream( 0 )
#endif
  !
  RETURN
  !
END SUBROUTINE close_files
