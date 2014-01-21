!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE offset_atom_wfc( nat, offset )
  !----------------------------------------------------------------------------
  !
  ! For each Hubbard atom, compute the index of the projector in the list of
  ! atomic wavefunctions. IMPORTANT NOTICE: if there is more than one state
  ! with the chosen value of Hubbard_l, the one selected for U calculation is
  ! the last state with the given l, but with strictly positive occupation, 
  !
  USE uspp_param,       ONLY : upf
  USE noncollin_module, ONLY : noncolin
  USE ions_base,        ONLY : ityp
  USE basis,            ONLY : natomwfc
  USE ldaU,             ONLY : Hubbard_l, Hubbard_U, Hubbard_alpha
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)  :: nat
  !
  INTEGER, INTENT(OUT) :: offset(nat)
  !
  INTEGER  :: counter, na, nt, n
  !
  !
  counter = 0
  offset(:) = -99
  !
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     DO n = 1, upf(nt)%nwfc
        !
        IF ( upf(nt)%oc(n) >= 0.D0 ) THEN
           !
           IF ( noncolin ) THEN
              !
              IF ( upf(nt)%has_so ) THEN
                 !

                IF (upf(nt)%oc(n)>0.D0.AND.upf(nt)%lchi(n)==Hubbard_l(nt).and.offset(na).eq.-99) &
                 offset(na) = counter

                 counter = counter + 2 * upf(nt)%lchi(n)
                 !
                 IF ( ABS( upf(nt)%jchi(n)-upf(nt)%lchi(n) - 0.5D0 ) < 1.D-6 ) &
                    counter = counter + 2
                 !
              ELSE
                 !
                IF (upf(nt)%oc(n)>0.D0.AND.upf(nt)%lchi(n)==Hubbard_l(nt)) &
                 offset(na) = counter
                 counter = counter + 2 * ( 2 * upf(nt)%lchi(n) + 1 )
                 !
              END IF
              !
           ELSE
              !
              IF ( upf(nt)%oc(n) > 0.D0 .AND. upf(nt)%lchi(n) == Hubbard_l(nt) ) &
                 offset(na) = counter
              !
              counter = counter + 2 * upf(nt)%lchi(n) + 1
              !
           END IF
        END IF
     END DO

     IF ( (Hubbard_U(nt).NE.0.D0 .OR. Hubbard_alpha(nt).NE.0.D0 ) .AND. &
         offset(na) < 0 ) CALL errore('offset_atom_wfc', 'wrong offset', na)

  END DO
  !
  IF ( counter.NE.natomwfc ) &
     CALL errore ('offset_atom_wfc', 'wrong number of wavefunctions', 1)
  !
  RETURN
  !
END SUBROUTINE offset_atom_wfc
!
