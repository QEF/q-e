! Copyright (C) 2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)


!----------------------------------------------------------------------
subroutine wannier_clean()
  !----------------------------------------------------------------------
  !    
  ! ... This routine deallocates all dynamically allocated arrays for wannier calc and closes all open files
  !
  USE wannier_new, only: wan_in, wan_pot, wannier_energy, wannier_occ, pp, coef
  USE io_files
  USE buffers
  USE basis,      ONLY : swfcatom
  USE ldaU,       ONLY : lda_plus_u
  USE fixed_occ, ONLY : one_atom_occupations
  
  IMPLICIT NONE
  LOGICAL :: opnd
  
  if(allocated(wan_in)) deallocate(wan_in)
  if(allocated(wan_pot)) deallocate(wan_pot)
  if(allocated(wannier_energy)) deallocate(wannier_energy)
  if(allocated(wannier_occ)) deallocate(wannier_occ)
  if(allocated(pp)) deallocate(pp)
  if(allocated(coef)) deallocate(coef)

  CALL close_buffer( iunwpp, 'keep' )
  CALL close_buffer( iunwf, 'keep' )
  
  IF ( .NOT. ( lda_plus_u .OR. one_atom_occupations ) ) THEN
     INQUIRE( UNIT = iunsat, OPENED = opnd )  
     IF ( opnd ) CALL close_buffer( iunsat, 'delete' )
  END IF
  
  IF(ALLOCATED(swfcatom)) DEALLOCATE(swfcatom)

  return
  !
end subroutine wannier_clean

