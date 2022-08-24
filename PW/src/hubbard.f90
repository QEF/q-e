!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------
SUBROUTINE hub_summary()
    !--------------------------------------------------------------------------
    !!
    !! This routine writes a summary about Hubbard parameters
    !!
    USE ldaU,               ONLY : lda_plus_u_kind, Hubbard_projectors, is_hubbard, &
                                   Hubbard_U, Hubbard_J0, Hubbard_U2, Hubbard_alpha, &
                                   Hubbard_beta, Hubbard_alpha_back, is_hubbard_back,  &
                                   Hubbard_J, Hubbard_l, Hubbard_n, lda_plus_u 
    USE noncollin_module,   ONLY : lspinorb
    USE ions_base,          ONLY : ntyp => nsp
    USE io_global,          ONLY : stdout
    !
    IMPLICIT NONE
    INTEGER :: nt
    !
    WRITE(stdout,'(5x,a)') 'Hubbard projectors: ' // TRIM(Hubbard_projectors)
    IF (lda_plus_u_kind == 0) THEN
       WRITE( stdout, '(5x,"Hubbard parameters of DFT+U (Dudarev formulation) in eV:")')
       DO nt = 1, ntyp
          IF (is_hubbard(nt)) THEN
            CALL write_hub_param (nt, Hubbard_U(nt), 'U', 1)
            CALL write_hub_param (nt, Hubbard_J0(nt), 'J0', 1)
            CALL write_hub_param (nt, Hubbard_alpha(nt), 'alpha', 1)
            CALL write_hub_param (nt, Hubbard_beta(nt), 'beta', 1)
          ENDIF
          IF (is_hubbard_back(nt)) THEN
            CALL write_hub_param (nt, Hubbard_U2(nt), 'U', 2)
            CALL write_hub_param (nt, Hubbard_alpha_back(nt), 'alpha', 2)
        ENDIF
       ENDDO
    ELSEIF (lda_plus_u_kind == 1) THEN
       WRITE( stdout, '(5x,"Hubbard parameters of DFT+U (Liechtenstein formulation) in eV:")')
       DO nt = 1, ntyp
          IF (Hubbard_U(nt) /= 0.d0) THEN
             IF (Hubbard_l(nt) == 0) THEN
                CALL write_hub_param (nt, Hubbard_U(nt), 'U', 1)
             ELSEIF (Hubbard_l(nt) == 1) THEN
                CALL write_hub_param (nt, Hubbard_U(nt), 'U', 1)
                CALL write_hub_param (nt, Hubbard_J(1,nt), 'J', 1)
             ELSEIF (Hubbard_l(nt) == 2) THEN
                CALL write_hub_param (nt, Hubbard_U(nt), 'U', 1)
                CALL write_hub_param (nt, Hubbard_J(1,nt), 'J', 1)
                CALL write_hub_param (nt, Hubbard_J(2,nt), 'B', 1)
             ELSEIF (Hubbard_l(nt) == 3) THEN
                CALL write_hub_param (nt, Hubbard_U(nt), 'U', 1)
                CALL write_hub_param (nt, Hubbard_J(1,nt), 'J', 1)
                CALL write_hub_param (nt, Hubbard_J(2,nt), 'E2', 1)
                CALL write_hub_param (nt, Hubbard_J(3,nt), 'E3', 1)
             ENDIF
          ENDIF
       ENDDO
       IF (lspinorb) THEN
           WRITE(stdout, '(5x,"DFT+U on averaged j=l+1/2,l-1/2 radial WFs")')
       ENDIF
    ELSEIF (lda_plus_u_kind == 2) THEN
       ! Info about the Hubbard V is printed by the routine alloc_neighborhood
       IF (ANY(Hubbard_J0(:)/=0.d0) .OR. ANY(Hubbard_alpha(:)/=0.d0) .OR. &
           ANY(Hubbard_beta(:)/=0.d0) .OR. ANY(Hubbard_alpha_back(:)/=0.d0)) &
           WRITE( stdout, '(5x,"Hubbard parameters of DFT+U+V (Dudarev formulation) in eV:")')
       DO nt = 1, ntyp
          IF (is_hubbard(nt)) THEN
            CALL write_hub_param (nt, Hubbard_J0(nt), 'J0', 1)
            CALL write_hub_param (nt, Hubbard_alpha(nt), 'alpha', 1)
            CALL write_hub_param (nt, Hubbard_beta(nt), 'beta', 1)
          ENDIF
          IF (is_hubbard_back(nt)) THEN
            CALL write_hub_param (nt, Hubbard_alpha_back(nt), 'alpha', 2)
        ENDIF
       ENDDO
    ENDIF
    !
    WRITE(stdout, '(/5x,"Internal variables: lda_plus_u =",l, ", lda_plus_u_kind = ",i1)') &
            lda_plus_u, lda_plus_u_kind
    !
    RETURN
    !
END SUBROUTINE hub_summary

!------------------------------------------------------------------------------
SUBROUTINE write_hub_param (nt, hub_parameter, hub_name, flag)
    !--------------------------------------------------------------------------
    !
    USE kinds,          ONLY : DP
    USE ions_base,      ONLY : atm
    USE constants,      ONLY : rytoev
    USE io_global,      ONLY : stdout
    USE ldaU,           ONLY : Hubbard_n, Hubbard_l, Hubbard_n2, Hubbard_l2, &
                               Hubbard_n3, Hubbard_l3, backall
    !
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nt
    CHARACTER(len=*), INTENT(IN) :: hub_name
    INTEGER, INTENT(IN) :: flag ! 1: first Hubbard channel
                                ! 2: second Hubbard channel
    REAL(DP) :: hub_parameter
    CHARACTER(LEN=1), EXTERNAL :: l_to_spdf
    !
    IF (hub_parameter /= 0.d0) THEN
       IF (flag==1) THEN
          WRITE(stdout,'(5x,a,i1,a,f8.4)')  &
              hub_name // '(' // TRIM(atm(nt)) // '-', Hubbard_n(nt), &
              l_to_spdf(Hubbard_l(nt),.FALSE.) // ') =', hub_parameter*rytoev
       ELSEIF (flag==2) THEN
          IF (.NOT.backall(nt)) THEN
             ! In this case there is one Hubbard channel for background states.
             WRITE(stdout,'(5x,a,i1,a,f8.4)')  &
                hub_name // '(' // TRIM(atm(nt)) // '-', Hubbard_n2(nt), &
                l_to_spdf(Hubbard_l2(nt),.FALSE.) // ') =', hub_parameter*rytoev
          ELSE
             ! In this case there are two Hubbard channels for background states.
             WRITE(stdout,'(5x,a,i1,a,i1,a,f8.4)')  &
                hub_name // '(' // TRIM(atm(nt))  // '-', Hubbard_n2(nt), &
                l_to_spdf(Hubbard_l2(nt),.FALSE.) // '-', Hubbard_n3(nt), &
                l_to_spdf(Hubbard_l3(nt),.FALSE.) // ') =', hub_parameter*rytoev
          ENDIF
       ENDIF
    ENDIF
    !
    RETURN
    !
END SUBROUTINE write_hub_param

!------------------------------------------------------------------------------
SUBROUTINE determine_hubbard_occ ( nt, lflag )
    !-----------------------------------------------------------------------
    !!
    !! This routine sets the starting Hubbard occupations based on the data
    !! read from the pseudopotentials.
    !! lflag=1 - first Hubbard manifold (main)
    !! lflag=2 - second Hubbard manifold
    !! lflag=3 - third Hubbard manifold
    !!
    USE kinds,        ONLY : DP
    USE uspp_param,   ONLY : upf
    USE ldaU,         ONLY : is_hubbard, is_hubbard_back, Hubbard_n, Hubbard_l, &
                             Hubbard_n2, Hubbard_l2, Hubbard_n3, Hubbard_l3,    &
                             hubbard_occ
    USE io_global,    ONLY : stdout
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nt          ! atomic type
    INTEGER, INTENT(IN) :: lflag       ! Hubbard channel
    !
    CHARACTER(LEN=2), ALLOCATABLE :: label(:)
    CHARACTER(LEN=2) :: label_aux
    CHARACTER(LEN=2) :: label_hub
    INTEGER :: i, & ! runs over all pseudo-atomic orbitals for the atomic type nt
               ldim
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    CHARACTER(LEN=1), EXTERNAL :: l_to_spdf
    CHARACTER(LEN=1), EXTERNAL :: capital
    LOGICAL :: first
    !
    IF ( upf(nt)%nwfc < 1 ) THEN
       CALL errore('determine_hubbard_occ', 'no atomic wavefunctions in &
               &pseudopotential file for species #' // upf(nt)%psd // new_line('a') // &
               &'use a pseudopotential file with atomic wavefunctions!', 1)
    ENDIF
    !
    IF (lflag==1) THEN
       label_hub = TRIM(int_to_char(Hubbard_n(nt))) // l_to_spdf(Hubbard_l(nt),.TRUE.)
    ELSEIF (lflag==2) THEN
       label_hub = TRIM(int_to_char(Hubbard_n2(nt))) // l_to_spdf(Hubbard_l2(nt),.TRUE.)
    ELSEIF (lflag==3) THEN
       label_hub = TRIM(int_to_char(Hubbard_n3(nt))) // l_to_spdf(Hubbard_l3(nt),.TRUE.)
    ELSE
       CALL errore('determine_hubbard_occ','Not allowed value of lflag',lflag)
    ENDIF
    !
    ldim = upf(nt)%nwfc
    ALLOCATE(label(ldim))
    !
    first=.true.
    DO i = 1, ldim
       ! Label of the i-th atomic orbital for the atomic type nt
       ! (if lowercase, then capitalize)
       label_aux = upf(nt)%els(i)
       label(i) = label_aux(1:1) // capital(label_aux(2:2))
       IF (label(i)==label_hub) THEN
          IF (first) THEN
             hubbard_occ(nt,lflag) = upf(nt)%oc(i)     
             first=.false.
          ELSE
             ! In FR-PPs the same Hubbard manifold appears twice,
             ! so we need to sum up all the electrons
             hubbard_occ(nt,lflag) = hubbard_occ(nt,lflag) &
                                    + upf(nt)%oc(i)
          ENDIF
       ENDIF
    ENDDO
    !
    IF (hubbard_occ(nt,lflag) < 0.d0) THEN
       WRITE(stdout,'(5x,a,8(1x,a))') TRIM(upf(nt)%psd) // &
               " pseudopotential contains the orbitals: ", (label(i),i=1,ldim)
       WRITE(stdout,'(5x,2a)') "Requested Hubbard manifold from the input: ", label_hub
       CALL errore('determine_hubbard_occ',&
               &'Mismatch between the requested and available manifolds',1)
    ENDIF
    !
    DEALLOCATE (label)
    !
    RETURN
    !
END SUBROUTINE determine_hubbard_occ
