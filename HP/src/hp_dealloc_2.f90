!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine hp_dealloc_2
  !----------------------------------------------------------------------
  !
  ! Deallocate various variables from the HP calculation
  !
  !USE ldaU,         ONLY : dist_s, ityp_s
  USE ldaU_hp,       ONLY : todo_atom, perturbed_atom, &
                            chi0, chi, ns, magn, ityp_new
  !
  IMPLICIT NONE
  !
  IF (ALLOCATED(todo_atom))       DEALLOCATE(todo_atom)
  IF (ALLOCATED(perturbed_atom))  DEALLOCATE(perturbed_atom)
  IF (ALLOCATED(chi0))            DEALLOCATE(chi0)
  IF (ALLOCATED(chi))             DEALLOCATE(chi)
  IF (ALLOCATED(ns))              DEALLOCATE(ns)
  IF (ALLOCATED(magn))            DEALLOCATE(magn)
  IF (ALLOCATED(ityp_new))        DEALLOCATE(ityp_new)
  !IF (ALLOCATED(dist_s))         DEALLOCATE(dist_s)
  !IF (ALLOCATED(ityp_s))         DEALLOCATE(ityp_s)
  !
  RETURN
  !
end subroutine hp_dealloc_2
