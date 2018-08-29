!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_find_inequiv_sites()
  !---------------------------------------------------------------------
  !
  ! Find inequivalent sites which will be perturbed (one at a time) 
  ! in the linear-response calculation.
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE ions_base,     ONLY : nat, tau, ityp, atm, ntyp => nsp
  USE symm_base,     ONLY : nsym, irt, s
  USE ldaU,          ONLY : is_hubbard
  USE ldaU_hp,       ONLY : at_equiv_criterium, todo_atom,             &
                            skip_atom, skip_type, merge_type,          &
                            perturb_only_atom, ns, nath_pert, atm_new, &
                            ityp_new, ntyp_new, disable_type_analysis
  !
  IMPLICIT NONE
  !
  INTEGER :: na, nb, nt, isym, i, counter
  ! counters on atoms and symmetry
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  ! Initialize labels for types of atoms
  !
  DO nt = 1, ntyp
     atm_new(nt) = atm(nt)
  ENDDO 
  !
  ! Allocate an auxiliary array for types of atoms
  !
  ALLOCATE(ityp_new(nat))
  ityp_new(:) = ityp(:) 
  !
  ! Determine which Hubbard atoms must be perturbed
  !
  IF ( at_equiv_criterium == 1 ) THEN
     CALL select_pert_based_on_occupations()
  ELSEIF ( at_equiv_criterium == 2 ) THEN
     CALL select_pert_based_on_type()
  ELSEIF ( at_equiv_criterium == 3 ) THEN
     CALL select_pert_based_on_sym()
  ELSE
     CALL errore ('hp_find_inequiv_sites', 'Not allowed value of at_equiv_criterium', 1)
  ENDIF
  !
  ! Check whether there is at least one Hubbard atom which must be perturbed
  !
  IF ( ALL(todo_atom(:) .eqv. .false.) ) CALL errore ('hp_find_inequiv_sites', &
             & 'There are no Hubbard atoms to perturb', 1) 
  !
  DO nt = 1, ntyp
     IF ( at_equiv_criterium.NE.1 .AND. skip_type(nt) .AND. merge_type(nt)==0 ) &
      & CALL errore ('hp_find_inequiv_sites', 'merge_type was not specified', 1)
  ENDDO
  !
  ! If the user requested to skip perturbing some specific atoms
  ! or atoms of some specific type, then inform the code about this.
  ! Warning: Make sure you know what you are doing! Use this option
  ! only when you know that the atom which you do not want to perturb
  ! is equivalent to some other atom (but the code does not recognize
  ! this from the symmetry analysis)
  !
  counter = 0
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     IF ( at_equiv_criterium.NE.1 .AND. todo_atom(na) &
          & .AND. (skip_atom(na).OR.skip_type(nt)) ) todo_atom(na) = .false.
     !
     ! Consider only one atom  
     !
     IF (perturb_only_atom(na)) THEN
        !
        todo_atom(:)  = .false.
        todo_atom(na) = .true. 
        !
        counter = counter + 1
        IF (counter > 1) THEN
           CALL errore('hp_find_inequiv_sites', &
             & 'More than one perturb_only_atom(na)=.true. not allowed',1)
        ENDIF
        !
     ENDIF
     !
  ENDDO
  !
  ! Count how many Hubbard atoms will be perturbed
  !
  nath_pert = 0
  DO na = 1, nat
     IF (todo_atom(na)) nath_pert = nath_pert + 1
  ENDDO
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE select_pert_based_on_occupations
  !
  ! Determine which Hubbard atoms must be perturbed by
  ! analyzing the unperturbed occupations.
  !
  USE ldaU_hp,        ONLY : docc_thr
  !
  IMPLICIT NONE
  LOGICAL, ALLOCATABLE :: done_type(:)
  INTEGER :: nt1, nt2, nb_start
  !
  IF (ANY(skip_type(:))) &
      & CALL errore ('hp_find_inequiv_sites', &
      & 'skip_type must not be setup from the input when at_equiv_criterium=1', 1)
  IF (ANY(merge_type(:).NE.0)) &
      & CALL errore ('hp_find_inequiv_sites', &
      & 'merge_type must not be setup from the input when at_equiv_criterium=1', 1)
  IF (ANY(skip_atom(:))) &
      & CALL errore ('hp_find_inequiv_sites', &
      & 'skip_atom cannot be used when at_equiv_criterium=1', 1)
  !
  ALLOCATE(done_type(ntyp))
  done_type(:) = .false.
  todo_atom(:) = .false.
  skip_atom(:) = .false.
  !
  ntyp_new = ntyp
  !
  IF ( nat == 1 ) THEN
     nt = ityp(1)
     IF (is_hubbard(nt)) todo_atom(1) = .true.
     RETURN
  ENDIF
  !
  DO na = 1, nat
     !
     nt1 = ityp(na)
     !
     IF ( is_hubbard(nt1) .AND. .NOT.skip_atom(na)) THEN
        !
        todo_atom(na)  = .true.
        !
        IF ( .NOT.disable_type_analysis ) THEN
           !
           ! If there are atoms of the same type but have different
           ! occupations, then keep track of this information and then
           ! use it at the post-procesing stage, which will be important
           ! when doing the averaging and reconstruction of the matrix
           ! elements of the response matrices.
           !
           IF (.NOT.done_type(nt1)) THEN
              done_type(nt1) = .true.
           ELSE
              ntyp_new = ntyp_new + 1
              ityp_new(na) = ntyp_new
              atm_new(ntyp_new) = TRIM(atm(nt1)) // TRIM(int_to_char(ntyp_new))
              nt1 = ntyp_new
           ENDIF
           !
        ENDIF
        !
        IF (na < nat) THEN
           nb_start = na+1
        ELSE
           nb_start = nat
        ENDIF
        !
        DO nb = nb_start, nat
           !
           nt2 = ityp(nb)
           !  
           IF ( is_hubbard(nt2) ) THEN
              !
              ! If there are atoms of different type but with equal
              ! occupations, then do not perturb these atoms and 
              ! keep track of which atom it corresponds to (by type).
              ! 
              IF ( ABS(ns(nb)-ns(na)).LT.docc_thr ) THEN
                 skip_atom(nb) = .true.
                 ityp_new(nb) = nt1
              ENDIF
              ! 
           ENDIF
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE select_pert_based_on_occupations

SUBROUTINE select_pert_based_on_type()
  !
  ! Pick up one representative Hubbard atom for each atomic type.
  ! Warning: atoms which have the same type but which are inequivalent
  ! by symmetry will not be distinguished in this case.
  !
  IMPLICIT NONE
  LOGICAL, ALLOCATABLE :: done_type(:)
  !
  todo_atom(:) = .false.
  !
  ALLOCATE(done_type(ntyp))
  done_type(:) = .false.
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     IF ( is_hubbard(nt) .AND. .NOT.done_type(nt) ) THEN
          !
          todo_atom(na) = .true.
          done_type(nt) = .true.
          !
     ENDIF
     !
  ENDDO
  !
  DEALLOCATE(done_type)
  !
  RETURN
  !
END SUBROUTINE select_pert_based_on_type

SUBROUTINE select_pert_based_on_sym()
  !
  ! Find which Hubbard atoms have the same type and are equivalent by symmetry.
  ! Find how many inequivalent Hubbard atoms there are in total 
  ! (i.e. how many perturbations to be considered).
  !
  IMPLICIT NONE
  !
  ! Start from the assumption that all Hubbard atoms must be perturbed
  !
  todo_atom(:) = .false.
  !
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) todo_atom(na) = .true.
  ENDDO
  !
  IF ( nat == 1 ) RETURN
  !
  DO na = 1, nat-1
     !
     DO nb = na+1, nat
        !
        IF ( todo_atom(na) .AND. todo_atom(nb) ) THEN
           !
           ! Check if the rotated atom nb coincides with the atom na, 
           ! and if they have the same type (see PW/src/checksym -> irt).
           !
           DO isym = 1, nsym
              IF ( irt(isym,nb) == na ) THEN
                 todo_atom(nb) = .false.
                 GO TO 10
              ENDIF
           ENDDO
           !
        ENDIF
        !
10 CONTINUE
        !
     ENDDO
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE select_pert_based_on_sym

END SUBROUTINE hp_find_inequiv_sites
