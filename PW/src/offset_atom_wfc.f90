!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE offset_atom_wfc( Hubbard_only, l_back, offset, counter )
  !----------------------------------------------------------------------------
  !
  ! For each Hubbard atom, compute the index of the projector in the list of
  ! atomic wavefunctions. IMPORTANT NOTICE: if there is more than one state
  ! with the chosen value of Hubbard_l, the one selected for U calculation is
  ! the last state with the given l, but with strictly positive occupation, 
  ! If Hubbard_only=.true., find offsets with respect to the array of atomic
  ! wavefunctions including only those with a Hubbard U term (wfcU in ldaU.f90)
  !
  USE uspp_param,       ONLY : upf
  USE noncollin_module, ONLY : noncolin
  USE ions_base,        ONLY : nat, ityp
  USE ldaU,             ONLY : Hubbard_l, Hubbard_l_back, Hubbard_l1_back, &
                               is_hubbard, is_hubbard_back, backall
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN)  :: Hubbard_only
  INTEGER, INTENT(OUT) :: offset(nat), counter
  INTEGER, INTENT(IN)  :: l_back
  !
  INTEGER  :: na, nt, n, l
  LOGICAL  :: hubbard_wfc, hubbard_wfc_b, hubbard_wfc_b1
  !
  !
  counter = 0
  offset(:) = -99
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     DO n = 1, upf(nt)%nwfc
        !
        hubbard_wfc    = .FALSE.
        hubbard_wfc_b  = .FALSE.
        hubbard_wfc_b1 = .FALSE.
        !
        IF ( upf(nt)%oc(n) >= 0.D0 ) THEN
           !
           l = upf(nt)%lchi(n)
           IF (is_hubbard(nt)) hubbard_wfc = ( upf(nt)%oc(n)>0.D0 .AND. l == Hubbard_l(nt) )
           IF (is_hubbard_back(nt)) THEN
              hubbard_wfc_b = ( upf(nt)%oc(n)>0.D0 .AND. l == Hubbard_l_back(nt) ) 
              IF (backall(nt)) THEN
                 hubbard_wfc_b1 = ( upf(nt)%oc(n)>0.D0 .AND. l == Hubbard_l1_back(nt) )
              ENDIF
           ENDIF
           !
           IF ( noncolin ) THEN
              !
              IF ( upf(nt)%has_so ) THEN
                 !
                 ! offset to be set at the first occurrence of required l
                 IF (hubbard_wfc .AND. offset(na).eq.-99 ) offset(na) = counter
                 !
                 IF (hubbard_wfc .OR. .NOT. hubbard_only) THEN
                    ! j = l-1/2, degeneracy 2l
                    counter = counter + 2*l
                    ! j = l+1/2, degeneracy 2*l+2
                    IF (ABS( upf(nt)%jchi(n)-l-0.5D0 ) < 1.D-6) &
                      counter = counter + 2
                 END IF
                 !
                 IF (hubbard_wfc .AND. hubbard_only) THEN
                    counter = counter + 2*l + 2
                 END IF
                 !
              ELSE
                 !
                 IF (hubbard_wfc) offset(na) = counter
                 !
                 IF (hubbard_wfc .OR. .NOT. hubbard_only) THEN
                    counter = counter + 2*( 2*l + 1 )
                 END IF
                 !
              END IF
              !
           ELSE
              !
              IF ( ( (hubbard_wfc   .AND. l_back.EQ.1) .OR. &
                     (hubbard_wfc_b .AND. l_back.EQ.2) .OR. &
                     (hubbard_wfc_b1.AND. l_back.EQ.3) )    &
                      .AND. offset(na).EQ.-99 ) offset(na) = counter
              !
              IF (hubbard_wfc .OR. hubbard_wfc_b .OR. hubbard_wfc_b1 .OR. &
                 .NOT. hubbard_only) THEN
                 counter = counter + 2*l + 1
              ENDIF
              !
           ENDIF
        ENDIF
     ENDDO
     !
     IF (is_hubbard(nt) .OR. is_hubbard_back(nt)) THEN
        !
        IF ( l_back.EQ.1 .AND. is_hubbard(nt) .AND. offset(na) < 0 ) &
            CALL errore('offset_atom_wfc', 'wrong offset', na)
        IF ( l_back.EQ.2 .AND. is_hubbard_back(nt) .AND. offset(na) < 0 ) &
            CALL errore('offset_atom_wfc', 'wrong offset back', na)
        IF ( l_back.EQ.3 .AND. backall(nt) .AND. offset(na) < 0 ) &
            CALL errore('offset_atom_wfc', 'wrong offset back1', na)
        !
     ENDIF      
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE offset_atom_wfc
!
