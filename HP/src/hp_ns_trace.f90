!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_ns_trace
  !-----------------------------------------------------------------------
  !
  ! This routine computes the trace of occupations and the magnetization.
  ! The unperturbed occupation matrices were read via hub_readin which
  ! calls read_file.
  !
  USE kinds,         ONLY : DP
  USE scf,           ONLY : rho
  USE lsda_mod,      ONLY : nspin
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp
  USE ldaU,          ONLY : Hubbard_l, is_hubbard, lda_plus_u_kind 
                            !ldim_u, nsg, neighood ! DFT+U+V
  USE ldaU_hp,       ONLY : ns, magn
  !
  IMPLICIT NONE
  INTEGER :: is, na, na1, na2, nt, nt1, nt2, m1, m2, ldim, viz
  REAL(DP), ALLOCATABLE :: nsaux(:,:) ! auxiliary array for occupations
  !
  ALLOCATE(ns(nat)) 
  ALLOCATE(nsaux(nat,nspin))
  ns(:)      = 0.0d0
  nsaux(:,:) = 0.0d0
  !
  IF (nspin==2) THEN
     ALLOCATE(magn(nat))
     magn(:) = 0.0d0
  ENDIF
  !
  IF (lda_plus_u_kind.EQ.0) THEN
     !   
     DO na = 1, nat
        !
        nt = ityp(na)
        !
        IF (is_hubbard(nt)) THEN
           !
           ldim = 2 * Hubbard_l(nt) + 1
           !ldim = ldim_u(nt)
           !
           DO is = 1, nspin
              DO m1 = 1, ldim
                 nsaux(na,is) = nsaux(na,is) + rho%ns(m1,m1,is,na)
              ENDDO
           ENDDO
           !
           IF (nspin==1) THEN
              ns(na) = 2.0d0 * nsaux(na,1) 
           ELSE
              ns(na)   = nsaux(na,1) + nsaux(na,2)
              magn(na) = nsaux(na,1) - nsaux(na,2)
           ENDIF
           !
        ENDIF
        !
     ENDDO
     !
  ELSEIF (lda_plus_u_kind.EQ.2) THEN
     !
     CALL errore("hp_ns_trace", 'The HP code does not support lda_plus_u_kind=2',1)
     !
     !DO na1 = 1, nat
     !   nt1 = ityp(na1)
     !   IF (is_hubbard(nt1)) THEN
     !      ldim = ldim_u(nt1)
     !      DO viz = 1, neighood(na1)%num_neigh
     !         na2 = neighood(na1)%neigh(viz)
     !         IF (na2==na1) THEN
     !            DO is = 1, nspin
     !               DO m1 = 1, ldim
     !                  nsaux(na1,is) = nsaux(na1,is) + DBLE(nsg(na1,m1,viz,m1,is))
     !               ENDDO
     !            ENDDO
     !            IF (nspin==1) THEN
     !               ns(na1) = 2.0d0 * nsaux(na1,1)
     !            ELSE
     !               ns(na1)   = nsaux(na1,1) + nsaux(na1,2)
     !               magn(na1) = nsaux(na1,1) - nsaux(na1,2)
     !            ENDIF
     !            GO TO 10
     !         ENDIF
     !      ENDDO 
     !   ENDIF
     !10  CONTINUE
     !ENDDO 
     ! 
  ELSE
    CALL errore("hp_ns_trace",'This lda_plus_u_kind is not supported',1)
  ENDIF
  !
  DEALLOCATE(nsaux)
  !
  RETURN
  !
END SUBROUTINE hp_ns_trace
