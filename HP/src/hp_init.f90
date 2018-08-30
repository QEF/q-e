!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine hp_init()
  !----------------------------------------------------------------------
  !
  ! Setup various variables for the HP calculation.
  !
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE io_global,     ONLY : stdout
  USE lsda_mod,      ONLY : nspin
  USE ldaU,          ONLY : Hubbard_lmax, is_hubbard
  USE ldaU_hp,       ONLY : perturbed_atom, chi0, chi, todo_atom,     &
                            nath, nqsh, nath_sc, code, nq1, nq2, nq3, &
                            determine_num_pert_only
  !
  IMPLICIT NONE
  INTEGER :: na, nt
  !
  ALLOCATE (todo_atom(nat))
  ALLOCATE (perturbed_atom(nat))
  !
  ! Determine the number of q-points in the q-mesh without symmetry
  !
  nqsh = nq1*nq2*nq3
  !
  ! Determine the number of Hubbard atoms
  !
  nath = 0
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) nath = nath + 1
  ENDDO
  !
  ! Check if all Hubbard atoms are listed first and then
  ! all other atoms follow in the ATOMIC_POSITIONS card.
  !
  IF ( nath < nat ) THEN
     DO na = 1, nat
        nt = ityp(na)
        IF ( na > nath .AND. is_hubbard(nt)) THEN
           WRITE( stdout, '(/5x,"WARNING! All Hubbard atoms must be listed first in the " &
                             & "ATOMIC_POSITIONS card of PWscf")')
           WRITE( stdout, '(5x,"Stopping...")')
           CALL hp_stop_smoothly(.FALSE.)
        ENDIF
     ENDDO
  ENDIF
  !
  ! Determine the total number of Hubbard atoms in the virtual supercell
  !
  nath_sc = nath * nqsh
  !
  ! Find inequivalent sites
  !
  CALL hp_find_inequiv_sites()
  !
  IF (.NOT.determine_num_pert_only) THEN
     !
     ! The first dimension of chi0 and chi runs over all possible 
     ! real+virtual atoms (nath_sc), whereas the second dimension runs over 
     ! real atoms in the primitive cell which can be perturbed (nat).
     !
     ALLOCATE (chi0(nath_sc, nat))
     ALLOCATE (chi(nath_sc, nat))
     chi0(:,:) = 0.0d0
     chi(:,:) = 0.0d0
     !
  ENDIF
  !
  RETURN
  !
end subroutine hp_init
