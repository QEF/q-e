!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE close_files(lflag)
  !----------------------------------------------------------------------------
  !
  ! ... Close all files and synchronize processes for a new scf calculation.
  !
  USE ldaU,          ONLY : lda_plus_u, U_projection
  USE control_flags, ONLY : twfcollect, io_level
  USE fixed_occ,     ONLY : one_atom_occupations
  USE io_files,      ONLY : prefix, iunwfc, iunsat, &
                            iunhub, iunefield, iunefieldm, iunefieldp, &
                            iunwfc_exx
  USE buffers,       ONLY : close_buffer
  USE mp_images,     ONLY : intra_image_comm
  USE mp,            ONLY : mp_barrier
  USE wannier_new,   ONLY : use_wannier
  USE bp,            ONLY : lelfield
  !
  IMPLICIT NONE
  !
  LOGICAL, intent(in) :: lflag
  !
  LOGICAL :: opnd
  !
  !  ... close buffer/file containing wavefunctions: discard if
  !  ... wavefunctions are written in xml format, save otherwise
  !
  IF ( lflag .AND. (twfcollect .OR. io_level < 0 )) THEN
     CALL close_buffer ( iunwfc, 'DELETE' )
  ELSE
     CALL close_buffer ( iunwfc, 'KEEP' )
  END IF
  !
  ! ... close files associated with the EXX calculation
  !
  INQUIRE( UNIT = iunwfc_exx, OPENED = opnd )
  IF ( opnd ) CALL close_buffer ( iunwfc_exx, 'DELETE' )
  !
  ! ... iunsat contains the (orthogonalized) atomic wfcs * S
  ! ... iunhub as above, only for wavefcts having an associated Hubbard U
  !
  IF ( lda_plus_u .AND. (U_projection.NE.'pseudo') ) THEN
     IF ( io_level < 0 ) THEN
        CALL close_buffer ( iunhub,'DELETE' )
     ELSE
        CALL close_buffer ( iunhub,'KEEP' )
     END IF
  END IF
  IF ( use_wannier .OR. one_atom_occupations ) THEN
     IF ( io_level < 0 ) THEN
        CALL close_buffer ( iunsat,'DELETE' )
     ELSE
        CALL close_buffer ( iunsat,'KEEP' )
     END IF
  END IF
  !
  ! ... close unit for electric field if needed
  !
  IF ( lelfield ) THEN
     !
     IF ( io_level < 0 ) THEN
        CALL close_buffer ( iunefield, 'DELETE' )
        CALL close_buffer ( iunefieldm,'DELETE' )
        CALL close_buffer ( iunefieldp,'DELETE' )
     ELSE
        CALL close_buffer ( iunefield, 'KEEP' )
        CALL close_buffer ( iunefieldm,'KEEP' )
        CALL close_buffer ( iunefieldp,'KEEP' )
     ENDIF
     !
  END IF
  !
  CALL mp_barrier( intra_image_comm )  
  !
  RETURN
  !
END SUBROUTINE close_files
