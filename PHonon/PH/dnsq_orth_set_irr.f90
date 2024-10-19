!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------------------
SUBROUTINE dnsq_orth_set_irr(irr, ipert0)
    !--------------------------------------------------------------------------------------
    !! Set lr_has_dnsorth and dnsorth for phonons in irreducible representation irr
    !--------------------------------------------------------------------------------------
    !
    USE uspp,      ONLY : okvan
    USE ions_base, ONLY : nat
    USE lsda_mod,  ONLY : nspin
    USE ldaU,      ONLY : lda_plus_u, Hubbard_lmax
    USE ldaU_ph,   ONLY : dnsorth
    USE modes,     ONLY : npert
    USE ldaU_lr,   ONLY : lr_has_dnsorth, lr_dnsorth
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: irr
    !! irreducible representation
    INTEGER, INTENT(IN) :: ipert0
    !! Offset in perturbation index
    !
    INTEGER :: ldim, ipert
    !
    IF (lda_plus_u .AND. okvan) THEN
        lr_has_dnsorth = .TRUE.
        !
        ldim = 2 * Hubbard_lmax + 1
        ALLOCATE(lr_dnsorth(ldim, ldim, nspin, nat, npert(irr)))
        DO ipert = 1, npert(irr)
            lr_dnsorth(:, :, :, :, ipert) = dnsorth(:, :, :, :, ipert + ipert0)
        ENDDO
    ENDIF
    !
 END SUBROUTINE dnsq_orth_set_irr
 !-----------------------------------------------------------------------------------------
 !
 !----------------------------------------------------------------------------------------
 SUBROUTINE deallocate_dnsorth()
 !----------------------------------------------------------------------------------------
    USE uspp,    ONLY : okvan
    USE ldaU,    ONLY : lda_plus_u
    USE ldaU_lr, ONLY : lr_has_dnsorth, lr_dnsorth
    IMPLICIT NONE
    lr_has_dnsorth = .FALSE.
    IF (lda_plus_u .AND. okvan) THEN
        DEALLOCATE(lr_dnsorth)
    ENDIF
 !----------------------------------------------------------------------------------------
 END SUBROUTINE deallocate_dnsorth
 !----------------------------------------------------------------------------------------
 