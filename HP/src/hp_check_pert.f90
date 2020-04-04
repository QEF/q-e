!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE hp_check_pert(na)
  !----------------------------------------------------------------------
  !
  ! Check if the given atom na must be perterbed or not.
  !
  USE ions_base,     ONLY : ityp
  USE lsda_mod,      ONLY : nspin
  USE ldaU_hp,       ONLY : nah_pert, todo_atom, perturbed_atom
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na ! the atom under consideration
  !
  INTEGER :: nt ! dummy index for atomic type
  !
  ! None of atoms are perterbed
  perturbed_atom(:) = .false.
  !
  ! Determine the type of a given atom na
  nt = ityp(na)
  !
  IF ( todo_atom(na) ) THEN
     !
     ! Only atom na will be perturbed
     perturbed_atom(na) = .true.
     !
     ! Keep track of the site number of the perturbed atom
     nah_pert = na 
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE hp_check_pert
