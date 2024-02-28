!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE offset_atom_wfc( hubbard_only, lflag, offset, counter )
  !----------------------------------------------------------------------------
  !!
  !! For each Hubbard atom, compute the index of the projector in the list of
  !! atomic wavefunctions.
  !! If hubbard_only=.true., find offsets with respect to the array of atomic
  !! wavefunctions including only those with a Hubbard U term (wfcU in ldaU.f90) 
  !! IMPORTANT NOTICE: only strictly positive occupation (i.e. occ > 0.d0) is 
  !! accepted for the Hubbard manifold.
  !!
  USE uspp_param,       ONLY : upf
  USE noncollin_module, ONLY : noncolin
  USE ions_base,        ONLY : nat, ityp, atm
  USE io_global,        ONLY : stdout
  USE ldaU,             ONLY : Hubbard_l, Hubbard_n, Hubbard_l2, Hubbard_n2,        &
                               Hubbard_l3, Hubbard_n3, is_hubbard, is_hubbard_back, &
                               backall, hubbard_occ, hubbard_projectors
  USE upf_utils,        ONLY : l_to_spdf, capital
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN)  :: hubbard_only
  INTEGER, INTENT(OUT) :: offset(nat), counter
  INTEGER, INTENT(IN)  :: lflag ! index to distinguish between Hubbard manifolds
  !
  INTEGER  :: na, nt, n, l, ldim
  LOGICAL  :: hubbard_wfc, hubbard_wfc2, hubbard_wfc3
  CHARACTER(LEN=2) :: s,             &
                      label_hub,     & ! Label of the first Hubbard manifold
                      label_hub2,    & ! Label of the second Hubbard manifold
                      label_hub3       ! Label of the third Hubbard manifold
  CHARACTER(LEN=2) :: label_aux
  CHARACTER(LEN=2), ALLOCATABLE :: label(:)
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  counter = 0
  offset(:) = -1
  !
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     ldim = upf(nt)%nwfc
     !
     WRITE(s,'(i2)') nt
     IF ( ( is_hubbard(nt) .OR. is_hubbard_back(nt) .OR. &
            Hubbard_projectors=="ortho-atomic"      .OR. &
            Hubbard_projectors=="norm-atomic" ) .AND. ldim < 1 ) THEN
        CALL errore('offset_atom_wfc', 'no atomic wavefunctions in &
               &pseudopotential file for species #' // s // new_line('a') // &
               &'use a pseudopotential file with atomic wavefunctions!', lflag)
     ENDIF
     !
     label_hub=''; label_hub2=''; label_hub3=''
     !
     IF (is_hubbard(nt)) &
           label_hub = TRIM(int_to_char(Hubbard_n(nt))) // &
                            l_to_spdf(Hubbard_l(nt))
     !
     IF (is_hubbard_back(nt)) THEN
           label_hub2 = TRIM(int_to_char(Hubbard_n2(nt))) // &
                             l_to_spdf(Hubbard_l2(nt))
        IF (backall(nt)) &
           label_hub3 = TRIM(int_to_char(Hubbard_n3(nt))) // &
                             l_to_spdf(Hubbard_l3(nt))
     ENDIF
     !
     ALLOCATE(label(ldim))
     !
     DO n = 1, ldim
        !
        hubbard_wfc  = .FALSE.
        hubbard_wfc2 = .FALSE.
        hubbard_wfc3 = .FALSE.
        !
        ! Label of the n-th atomic orbital for the atom na of type nt
        ! (if lowercase, then capitalize)
        label_aux = upf(nt)%els(n)
        label(n) = label_aux(1:1) // capital(label_aux(2:2))
        !
        IF (TRIM(label(n))=='') &
                CALL errore('offset_atom_wfc', 'The pseudo for ' // atm(nt) // &
                & ' does not contain labels for atomic orbitals! &
                &Please add them by hand in the pseudo.',1)
        !
        IF ( upf(nt)%oc(n) >= 0.D0 ) THEN
           !
           l = upf(nt)%lchi(n)
           !
           ! Note: For FR-PPs, the j = l-1/2 is full, the j = l+1/2 is empty.
           ! For example, for Eu-4f one 4f orbital has non-zero occupancy, while 
           ! other 4f orbital is empty. However, hubbard_occ is defined as a sum
           ! of both occupied and empty 4f orbitals (see determine_hubbard_occ).
           ! Therefore, both occupied and empty 4f states will be considered in 
           ! the calculation (the full 4f shell), which is consistent with the 
           ! NR/SR case where the full 4f shell is considered. 
           !
           IF (is_hubbard(nt)) THEN
              IF (label(n)==label_hub) THEN
                 IF (hubbard_occ(nt,1) > 0.D0) THEN
                    hubbard_wfc = .TRUE.
                 ELSE
                    CALL errore('offset_atom_wfc', 'Hubbard manifold with &
                            &zero occupations is not allowed',1) 
                 ENDIF
              ENDIF
           ENDIF
           !
           IF (is_hubbard_back(nt)) THEN
              IF (label(n)==label_hub2) THEN
                 IF (hubbard_occ(nt,2) > 0.D0) THEN
                    hubbard_wfc2 = .TRUE.
                 ELSE
                    CALL errore('offset_atom_wfc', 'Hubbard manifold with &
                            &zero occupations is not allowed',1)
                 ENDIF
              ENDIF
              IF (backall(nt) .AND. (label(n)==label_hub3)) THEN
                 IF (hubbard_occ(nt,3) > 0.D0) THEN
                    hubbard_wfc3 = .TRUE.
                 ELSE
                    CALL errore('offset_atom_wfc', 'Hubbard manifold with &
                            &zero occupations is not allowed',1)
                 ENDIF
              ENDIF
           ENDIF
           !
           IF ( noncolin ) THEN
              !
              IF ( upf(nt)%has_so ) THEN
                 !
                 ! offset to be set at the first occurrence of required manifold
                 IF (hubbard_wfc .AND. offset(na).eq.-1) offset(na) = counter
                 !
                 IF (hubbard_wfc .OR. .NOT.hubbard_only) THEN
                    ! j = l-1/2, degeneracy 2l
                    counter = counter + 2*l
                    ! j = l+1/2, degeneracy 2*l+2
                    IF (ABS( upf(nt)%jchi(n)-l-0.5D0 ) < 1.D-6) &
                      counter = counter + 2
                 ENDIF
                 !
                 IF (hubbard_wfc .AND. hubbard_only) &
                    counter = counter + 2*l + 2
                 !
              ELSE
                 !
                 IF (hubbard_wfc) offset(na) = counter
                 !
                 IF (hubbard_wfc .OR. .NOT.hubbard_only) &
                    counter = counter + 2*( 2*l + 1 )
                 !
              ENDIF
              !
           ELSE
              !
              IF ( ( (hubbard_wfc  .AND. lflag.EQ.1) .OR. &
                     (hubbard_wfc2 .AND. lflag.EQ.2) .OR. &
                     (hubbard_wfc3 .AND. lflag.EQ.3) ) )  &
                      offset(na) = counter
              !
              IF (hubbard_wfc .OR. hubbard_wfc2 .OR. hubbard_wfc3 .OR. &
                 .NOT. hubbard_only) THEN
                 counter = counter + 2*l + 1
              ENDIF
              !
           ENDIF
        ENDIF
     ENDDO
     !
     IF ( ( is_hubbard(nt) .OR. is_hubbard_back(nt) ) .AND. offset(na) < 0 ) THEN
        WRITE(stdout,'(5x,a,8(1x,a))') TRIM(upf(nt)%psd) // &
                " pseudopotential contains the orbitals: ", (label(n),n=1,ldim)
        IF (lflag==1 .AND. is_hubbard(nt)) THEN
           WRITE(stdout,'(5x,2a)') "Requested Hubbard manifold from the input: ", label_hub
        ELSEIF (lflag==2 .AND. is_hubbard_back(nt)) THEN
           WRITE(stdout,'(5x,2a)') "Requested Hubbard manifold (2nd) from the input: ", label_hub2
        ELSEIF (lflag==3 .AND. is_hubbard_back(nt) .AND. backall(nt)) THEN
           WRITE(stdout,'(5x,2a)') "Requested Hubbard manifold (3rd) from the input: ", label_hub3
        ELSE
           GO TO 11 
        ENDIF
        CALL errore('offset_atom_wfc','Mismatch between the requested and available manifolds',lflag)
11      CONTINUE    
     ENDIF      
     !
     DEALLOCATE(label)
     !
  ENDDO
  !
  RETURN
  !
END SUBROUTINE offset_atom_wfc
!
