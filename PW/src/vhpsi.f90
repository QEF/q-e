!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine vhpsi (ldap, np, mps, psip, hpsi)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the Hubbard potential applied to the electronic
  ! structure of the current k-point. The result is added to hpsi
  !
  USE kinds,     ONLY : DP
  USE becmod,    ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
  USE ldaU,      ONLY : Hubbard_lmax, Hubbard_l, is_Hubbard, nwfcU, wfcU, &
                        offsetU
  USE lsda_mod,  ONLY : current_spin
  USE scf,       ONLY : v
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE control_flags, ONLY : gamma_only
  USE mp,        ONLY: mp_sum
  !
  implicit none
  !
  integer, intent (in) :: ldap, np, mps
  complex(DP), intent(in) :: psip (ldap, mps)
  complex(DP), intent(inout) :: hpsi (ldap, mps)
  !
  integer :: na, nt, ldim
  REAL(DP), ALLOCATABLE :: rtemp(:,:)
  COMPLEX(DP), ALLOCATABLE :: ctemp(:,:), vtemp(:,:)
  type (bec_type) :: proj

  CALL start_clock('vhpsi')
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
  !
  CALL allocate_bec_type ( nwfcU,mps, proj )
  CALL calbec (np, wfcU, psip, proj)
  !
  DO nt = 1, ntyp
     IF ( is_hubbard(nt) ) THEN
        ldim = 2*Hubbard_l(nt) + 1
        IF  (gamma_only) THEN
           ALLOCATE ( rtemp(ldim,mps) )
        ELSE
           ALLOCATE ( ctemp(ldim,mps), vtemp(ldim,ldim) )
        END IF
        DO na = 1, nat  
           IF ( nt == ityp (na) ) THEN
              IF (gamma_only) THEN
                 CALL DGEMM ('n','n', ldim,mps,ldim, 1.0_dp, &
                      v%ns(1,1,current_spin,na),2*Hubbard_lmax+1, &
                      proj%r(offsetU(na)+1,1),nwfcU, 0.0_dp, rtemp, ldim)
                 CALL DGEMM ('n','n', 2*np, mps, ldim, 1.0_dp, &
                     wfcU(1,offsetU(na)+1), 2*ldap, rtemp, ldim, &
                     1.0_dp, hpsi, 2*ldap)
              ELSE
                 vtemp(:,:) =  v%ns(1:ldim,1:ldim,current_spin,na)
                 CALL ZGEMM ('n','n', ldim,mps,ldim, (1.0_dp,0.0_dp), &
                      vtemp,ldim, proj%k(offsetU(na)+1,1),nwfcU, &
                      (0.0_dp, 0.0_dp), ctemp, ldim)
                 CALL ZGEMM ('n','n', np, mps, ldim, (1.0_dp,0.0_dp), &
                      wfcU(1,offsetU(na)+1), ldap, ctemp, ldim, &
                      (1.0_dp,0.0_dp), hpsi, ldap)
              ENDIF
           ENDIF
        ENDDO
        IF (gamma_only) THEN
           DEALLOCATE ( rtemp )
        ELSE
           DEALLOCATE ( vtemp, ctemp )
        ENDIF
     ENDIF
  ENDDO
  !
  CALL deallocate_bec_type (proj)
  !
  CALL stop_clock('vhpsi')

  RETURN

END subroutine vhpsi

subroutine vhpsi_nc (ldap, np, mps, psip, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Noncollinear version (A. Smogunov). 
  !
  USE kinds,            ONLY : DP
  USE ldaU,             ONLY : Hubbard_lmax, Hubbard_l, is_Hubbard, nwfcU, &
                               wfcU, offsetU
  USE scf,              ONLY : v
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp
  USE noncollin_module, ONLY : npol
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  !
  implicit none
  !
  integer, intent (in) :: ldap, np, mps
  complex(DP), intent(in) :: psip (ldap*npol, mps)
  complex(DP), intent(inout) :: hpsi (ldap*npol, mps)
  !
  integer :: ibnd, na, nwfc, is1, is2, nt, m1, m2
  complex(DP) :: temp, zdotc 
  complex(DP), allocatable :: proj(:,:)

  CALL start_clock('vhpsi')
  ALLOCATE( proj(nwfcU, mps) )

!-- FIXME: to be replaced with ZGEMM
! calculate <psi_at | phi_k> 
  DO ibnd = 1, mps
    DO na = 1, nwfcU
      proj(na, ibnd) = zdotc (ldap*npol, wfcU(1, na), 1, psip(1, ibnd), 1)
    ENDDO
  ENDDO
#if defined(__MPI)
  CALL mp_sum ( proj, intra_bgrp_comm )
#endif
!--

  do ibnd = 1, mps  
    do na = 1, nat  
       nt = ityp (na)  
       if ( is_hubbard(nt) ) then  
          nwfc = 2 * Hubbard_l(nt) + 1

          do is1 = 1, npol
           do m1 = 1, nwfc 
             temp = 0.d0
             do is2 = 1, npol
              do m2 = 1, nwfc  
                temp = temp + v%ns_nc( m1, m2, npol*(is1-1)+is2, na) * &
                              proj(offsetU(na)+(is2-1)*nwfc+m2, ibnd)
              enddo
             enddo
             call zaxpy (ldap*npol, temp, wfcU(1,offsetU(na)+(is1-1)*nwfc+m1),&
                         1, hpsi(1,ibnd),1)
           enddo
          enddo

       endif
    enddo
  enddo

  deallocate (proj)
  CALL stop_clock('vhpsi')

  return
end subroutine vhpsi_nc

