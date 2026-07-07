!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE dnsq_store(npe, imode0)
!----------------------------------------------------------------------------
  !! Store the computed dnsscf in the full matrix dnsscf_all_modes
  !! (i.e for all modes and not only for the npe irreducible representations)
  !
  USE ions_base,     ONLY : nat, ityp
  USE lsda_mod,      ONLY : nspin
  USE ldaU,          ONLY : is_hubbard, Hubbard_l
  USE ldaU_lr,       ONLY : dnsscf
  USE ldaU_ph,       ONLY : dnsscf_all_modes
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN) :: npe
  !! the number of perturbations
  INTEGER , INTENT(IN) :: imode0
  !! the position of the modes
  !
  INTEGER :: ipert, nah, nt, is, m1, m2
  !
  DO ipert = 1, npe
   DO nah = 1, nat
      nt = ityp(nah)
      IF (is_hubbard(nt)) THEN
         DO is = 1, nspin
            DO m1 = 1, 2*Hubbard_l(nt)+1
               DO m2 = 1, 2*Hubbard_l(nt)+1
                  dnsscf_all_modes(m1,m2,is,nah,imode0+ipert) = &
                                   dnsscf(m1,m2,is,nah,ipert)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDDO
END SUBROUTINE dnsq_store
!----------------------------------------------------------------------------
