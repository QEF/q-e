!
! Copyright (C) 2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE write_modes_out(irr, imode0)
!
! This routine writes the displacements on the representation irr that
! starts at mode imode0
!
USE kinds, ONLY : DP
USE ions_base, ONLY : nat
USE modes, ONLY : u, npert
USE io_global, ONLY : stdout
IMPLICIT NONE
INTEGER, INTENT(IN) :: imode0, irr
INTEGER :: mu, nu

WRITE( stdout, '(5x,"Irreps are as follows:",/)')
IF (npert (irr) .eq.1) THEN
   WRITE( stdout, '(20x," mode # ",i3)') imode0 + 1
   WRITE( stdout, '(20x," (",2f10.5,"   ) ")')  ( (u (mu, nu) ,&
         &nu = imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
ELSEIF (npert (irr) .eq.2) THEN
   WRITE( stdout, '(2(10x," mode # ",i3,16x))') imode0 + 1, &
         imode0 + 2
   WRITE( stdout, '(2(10x," (",2f10.5,"   ) "))')  ( (u (mu, nu) , nu &
         &= imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
ELSEIF (npert (irr) .eq.3) THEN
   WRITE( stdout, '(4x,3(" mode # ",i3,13x))') imode0 + 1, imode0 &
         + 2, imode0 + 3
   WRITE( stdout, '((5x,3("(",2f10.5," ) ")))') ( (u (mu, nu) , &
         nu = imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
ELSE
   WRITE( stdout, '(4x,4(" mode # ",i3,13x))') imode0 + 1, imode0 &
         + 2, imode0 + 4
   WRITE( stdout, '((5x,4("(",2f10.5," ) ")))') ( (u (mu, nu) , &
         nu = imode0 + 1, imode0 + npert (irr) ) , mu = 1, 3 * nat)
ENDIF

RETURN
END SUBROUTINE write_modes_out
