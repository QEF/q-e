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
  USE ldaU,          ONLY : lda_plus_u
  USE control_flags, ONLY : twfcollect, io_level
  USE io_files,      ONLY : prefix, iunwfc, iunigk, iunat, iunsat
  USE buffers,       ONLY : close_buffer
  USE mp_global,     ONLY : intra_image_comm
  USE mp,            ONLY : mp_barrier
  !
  IMPLICIT NONE
  !
  LOGICAL :: opnd
  !  ... close buffer/file containing wavefunctions: discard if
  !  ... wavefunctions are written in xml format, save otherwise
  !
  IF ( twfcollect .OR. io_level < 0 ) THEN
     CALL close_buffer ( iunwfc, 'DELETE' )
  ELSE
     CALL close_buffer ( iunwfc, 'KEEP' )
  END IF
  !
  ! ... iunigk is kept open during the execution - close and remove
  !
  INQUIRE( UNIT = iunigk, OPENED = opnd )
  IF ( opnd ) CLOSE( UNIT = iunigk, STATUS = 'DELETE' )
  !
  ! ... iunat  contains the (orthogonalized) atomic wfcs 
  ! ... iunsat contains the (orthogonalized) atomic wfcs * S
  !
  IF ( lda_plus_u ) THEN
     !
     INQUIRE( UNIT = iunat, OPENED = opnd )  
     IF ( opnd ) CLOSE( UNIT = iunat, STATUS = 'KEEP' )
     INQUIRE( UNIT = iunsat, OPENED = opnd )  
     IF ( opnd ) CLOSE( UNIT = iunsat, STATUS = 'KEEP' )
     !
  END IF
  !
  CALL mp_barrier( intra_image_comm )  
  !
#if defined (__T3E)
  !
  ! ... set streambuffers off
  !
  CALL set_d_stream( 0 )
#endif
  !
  RETURN
  !
END SUBROUTINE close_files
