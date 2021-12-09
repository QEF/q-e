!
! Copyright (C) 2002-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE engine_to_path_pos(idx)
  !-----------------------------------------------------------------------------
  !
  ! ... Copy atomic positions (tau) and atom types (ityp) read from file
  ! ... to array pos(:,idx) and typ(:), where idx is the index of image 
  ! ... Translate positions by lattice vectors to make the path smooth,
  ! ... verify that typ is the same array for all images
  !
  USE kinds,         ONLY : DP
  !
  USE path_input_parameters_module, ONLY : input_images, minimum_image
  USE path_input_parameters_module, ONLY : nat, alat 
  USE path_input_parameters_module, ONLY : pos, typ
  USE path_io_units_module,         ONLY : iunpath
  !
  USE ions_base, ONLY : tau, ityp
  USE cell_base, ONLY : bg, at
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: idx
  !
  INTEGER :: iat
  REAL(DP), ALLOCATABLE :: pos0(:,:), pos1(:,:)
  ! atomic positions (in crystal units) of the previous and current image
  !
  ! tau is already in the internal QE units
  !
  ! is this really necessary? (GS)
  if(.not.allocated(pos)) allocate(pos(3*nat,input_images))
  pos(:,idx) = 0.0_dp
  !
  ! ... note that this positions array is in Bohr
  !
  pos(1:3*nat,idx) = reshape( tau, (/ 3 * nat /) ) * alat
  !
  ! If requested, use the translational periodicity of the unit cell to ensure
  ! that the path is smooth (even if atoms in the input images do not move on
  ! a smooth path). It adopts a "minimum image criterion", which assumes that 
  ! atoms do not be displaced by more than half the unit cell size from one 
  ! input image to the next. If this happens, a periodic replica of that atom
  ! is chosen to avoid atomic "jumps". If not requested, just give a warning 
  ! when an atom moves that much. -GS
  !
  ALLOCATE( pos0(3,nat), pos1(3,nat) )
  !
  ! atomic positions in current image
  pos1 = reshape( pos(:,idx),   (/ 3, nat /) ) / alat
  CALL cryst_to_cart( nat, pos1(1,1), bg, -1 )
  ! refold them within the unit cell around the origin
  IF ( minimum_image ) pos1 = pos1(:,:) - anint(pos1(:,:))
  !
  IF ( idx > 1 ) THEN
     !
     ! atomic positions in previous image
     pos0 = reshape( pos(:,idx-1), (/ 3, nat /) ) / alat
     CALL cryst_to_cart( nat, pos0(1,1), bg, -1 )
     !
     DO iat = 1,nat
        !
        ! translate atom by a lattice vector if needed
        ! N.B.: this solves the problem only when |p1-p0|<1.0
        !
        IF ( minimum_image ) THEN
           WHERE( (pos1(:,iat) - pos0(:,iat)) > 0.5_DP )
              pos1(:,iat) = pos1(:,iat) - 1.0_DP
           ENDWHERE
           !
           WHERE( (pos1(:,iat) - pos0(:,iat)) < -0.5_DP ) 
              pos1(:,iat) = pos1(:,iat) + 1.0_DP
           ENDWHERE
        ELSE
           IF ( ANY(ABS(pos1(:,iat) - pos0(:,iat)) > 0.5_DP) ) THEN
              WRITE ( iunpath, '(/,5x,A,I5,A,I3,A,I3,/,5x,A)' ) "WARNING: atom", iat, &
                 " moved more than 1/2 alat from image", idx-1, " to image", idx, &
                 "You can set minimum_image to true to avoid jumps in the path"
           ENDIF
        ENDIF
     ENDDO
  ENDIF
  !
  IF ( minimum_image ) THEN
     !
     ! update positions only if requested
     CALL cryst_to_cart( nat, pos1(1,1), at, 1 )
     pos(1:3*nat,idx) = reshape( pos1, (/ 3 * nat /) ) * alat
     !
  ENDIF
  !
  DEALLOCATE( pos0, pos1 )
  !
  ! consistency check on atomic type, just to be sure... (GS)
  !
  IF ( idx == 1 ) THEN
     ! is this really necessary? (GS)
     if(.not.allocated(typ)) allocate(typ(nat))
     typ = ityp(:)
  ELSE
     IF ( ANY( typ .NE. ityp ) ) CALL errore("engine_to_path_pos", &
        "inconsistency of atomic species", idx )
  ENDIF
  !
  RETURN
  !
END SUBROUTINE engine_to_path_pos
!
