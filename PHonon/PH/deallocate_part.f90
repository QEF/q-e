!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------
subroutine deallocate_part()
!----------===============-------------------------

  USE partial, ONLY : comp_irr, done_irr, atomo
  USE el_phon, ONLY : done_elph, comp_elph
  IMPLICIT NONE

  if (allocated(comp_irr)) deallocate (comp_irr)
  if (allocated(done_irr)) deallocate (done_irr)
  if (allocated(comp_elph)) deallocate (comp_elph)
  if (allocated(done_elph)) deallocate (done_elph)
  if (allocated(atomo)) deallocate (atomo)

  return
end subroutine deallocate_part
