!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE close_files( lflag )
  !----------------------------------------------------------------------------
  !! Close all files and synchronize processes for a new scf calculation.
  !
  USE ldaU,          ONLY: lda_plus_u, Hubbard_projectors
  USE control_flags, ONLY: io_level
  USE fixed_occ,     ONLY: one_atom_occupations
  USE io_files,      ONLY: prefix, iunwfc, iunsat, &
                           iunhub, iunefield, iunefieldm, iunefieldp, &
                           iunwfc_exx, iunhub_noS
  USE buffers,       ONLY: close_buffer
  USE mp_images,     ONLY: intra_image_comm
  USE mp,            ONLY: mp_barrier
  USE wannier_new,   ONLY: use_wannier
  USE bp,            ONLY: lelfield
#if defined (__OSCDFT)
  USE plugin_flags,     ONLY : use_oscdft
  USE oscdft_base,      ONLY : oscdft_ctx
  USE oscdft_functions, ONLY : oscdft_close_files
#endif
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lflag
  !
  CHARACTER(LEN=6) :: close_option
  LOGICAL :: opnd
  !
  !  ... delete buffers 
  !  ... 1) at convergence, unless high disk I/O is required
  !  ... 2) always, if minimal disk I/O is required
  !
  IF ( (lflag .AND. io_level == 0) .OR. (io_level < 0) ) THEN
     close_option = 'DELETE'
  ELSE
     close_option = 'KEEP'
  END IF
  !
  !  ... close buffer/file containing wavefunctions
  !
  CALL close_buffer ( iunwfc, close_option )
  !
  ! ... close files associated with the EXX calculation
  !
  INQUIRE( UNIT = iunwfc_exx, OPENED = opnd )
  IF ( opnd ) CALL close_buffer( iunwfc_exx, 'DELETE' )
  !
  ! ... iunsat contains the (orthogonalized) atomic wfcs * S
  ! ... iunhub  as above, only for wfcs * S having an associated Hubbard U
  !
  IF ( lda_plus_u .AND. (Hubbard_projectors /= 'pseudo') ) THEN
     CALL close_buffer( iunhub, 'DELETE' )
     INQUIRE( UNIT = iunhub_noS, OPENED = opnd )
     IF ( opnd ) CALL close_buffer( iunhub_noS, close_option )
  END IF
  !
  IF ( use_wannier .OR. one_atom_occupations ) THEN
     CALL close_buffer( iunsat, close_option )
  END IF
  !
  ! ... close unit for electric field if needed
  !
  IF ( lelfield ) THEN
     CALL close_buffer( iunefield,  close_option )
     CALL close_buffer( iunefieldm, close_option )
     CALL close_buffer( iunefieldp, close_option )
  END IF
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) CALL oscdft_close_files(oscdft_ctx)
#endif
  !
  CALL mp_barrier( intra_image_comm )  
  !
  RETURN
  !
END SUBROUTINE close_files
