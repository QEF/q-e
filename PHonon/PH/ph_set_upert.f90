!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------------------------------------------------------------------
SUBROUTINE ph_set_upert_phonon(irr)
   !--------------------------------------------------------------------------------------
   !! Set lr_npert, upert, and upert_mp for phonons in irreducible representation irr
   !--------------------------------------------------------------------------------------
   !
   USE modes,        ONLY : npert, t, tmq
   USE lr_symm_base, ONLY : nsymq, minus_q, lr_npert, upert, upert_mq
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN) :: irr
   !! irreducible representation
   !
   INTEGER :: ipert, jpert
   !! Counter on perturbations
   INTEGER :: isym
   !! Counter on symmetries
   !
   ! Set symmetry representation in lr_symm_base
   !
   lr_npert = npert(irr)
   !
   ALLOCATE(upert(lr_npert, lr_npert, nsymq))
   !
   DO isym = 1, nsymq
      DO ipert = 1, lr_npert
         DO jpert = 1, lr_npert
            upert(jpert, ipert, isym) = t(jpert, ipert, isym, irr)
         ENDDO
      ENDDO
   ENDDO
   !
   IF (minus_q) THEN
      ALLOCATE(upert_mq(lr_npert, lr_npert))
      DO ipert = 1, lr_npert
         DO jpert = 1, lr_npert
            upert_mq(jpert, ipert) = tmq(jpert, ipert, irr)
         ENDDO
      ENDDO
   ENDIF ! minus_q
   !
END SUBROUTINE ph_set_upert_phonon
!-----------------------------------------------------------------------------------------
!
!-----------------------------------------------------------------------------------------
SUBROUTINE ph_set_upert_e()
   !--------------------------------------------------------------------------------------
   !! Set lr_npert, upert, and upert_mp for electric field perturbation.
   !--------------------------------------------------------------------------------------
   !
   USE symm_base,    ONLY : s
   USE lr_symm_base, ONLY : nsymq, minus_q, lr_npert, upert, upert_mq
   !
   IMPLICIT NONE
   !
   INTEGER :: ipol, jpol
   !! Counter on perturbations
   INTEGER :: isym
   !! Counter on symmetries
   !
   ! Set symmetry representation in lr_symm_base
   !
   lr_npert = 3
   !
   ALLOCATE(upert(lr_npert, lr_npert, nsymq))
   !
   DO isym = 1, nsymq
      DO ipol = 1, 3
         DO jpol = 1, 3
            upert(jpol, ipol, isym) = s(ipol, jpol, isym)
         ENDDO
      ENDDO
   ENDDO
   !
   IF (minus_q) THEN
      !
      ! upert_mq is the rotation matrix for symmetry S such that T * S * q = q + G.
      ! E field perturbation is applied only for q = 0, where T*q = q, i.e., S = identity.
      ! Thus, upert_mq is the 3*3 identity matrix.
      !
      ALLOCATE(upert_mq(lr_npert, lr_npert))
      upert_mq = (0.d0, 0.d0)
      DO ipol = 1, lr_npert
         upert_mq(ipol, ipol) = (1.d0, 0.d0)
      ENDDO
   ENDIF ! minus_q
   !
END SUBROUTINE ph_set_upert_e
!----------------------------------------------------------------------------------------
!
!----------------------------------------------------------------------------------------
SUBROUTINE ph_deallocate_upert()
!----------------------------------------------------------------------------------------
   USE lr_symm_base, ONLY : nsymq, minus_q, lr_npert, upert, upert_mq
   IMPLICIT NONE
   DEALLOCATE(upert)
   IF (minus_q) DEALLOCATE(upert_mq)
!----------------------------------------------------------------------------------------
END SUBROUTINE ph_deallocate_upert
!----------------------------------------------------------------------------------------
