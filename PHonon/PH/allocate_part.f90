!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine allocate_part ( nat )
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays for the control of partial computation
  ! of the dynamical matrix
  !
  USE partial, ONLY : comp_irr, done_irr, atomo
  USE el_phon, ONLY : comp_elph, done_elph, elph
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nat
  !
  !  allocate space for several arrays which control the run
  !
  allocate (comp_irr (  0:3 * nat))
  allocate (done_irr (  0:3 * nat))
  IF (elph) THEN
     allocate (comp_elph (1:3 * nat))
     allocate (done_elph (1:3 * nat))
  ENDIF
  allocate (atomo    (  nat))
  atomo(:) = 0
  return
end subroutine allocate_part
