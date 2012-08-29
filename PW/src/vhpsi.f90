!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
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
  ! of the current k-point, the result is added to hpsi
  !
  USE kinds,     ONLY : DP
  USE becmod,    ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
  USE ldaU,      ONLY : Hubbard_lmax, Hubbard_l, HUbbard_U, Hubbard_alpha, &
                        swfcatom, oatwfc
  USE lsda_mod,  ONLY : current_spin
  USE scf,       ONLY : v
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE basis,     ONLY : natomwfc
  USE control_flags, ONLY : gamma_only
  USE mp,        ONLY: mp_sum
  !
  implicit none
  !
  integer, intent (in) :: ldap, np, mps
  complex(DP), intent(in) :: psip (ldap, mps)
  complex(DP), intent(inout) :: hpsi (ldap, mps)
  !
  integer :: ibnd, na, nt, m1, m2
  complex(DP) :: temp
  type (bec_type) :: proj

  CALL start_clock('vhpsi')
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in oatwfc
  !
  call allocate_bec_type ( natomwfc,mps, proj )
  CALL calbec (np, swfcatom, psip, proj)
  do ibnd = 1, mps  
     do na = 1, nat  
        nt = ityp (na)  
        if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then  
           do m1 = 1, 2 * Hubbard_l(nt) + 1 
              temp = 0.d0
              if (gamma_only) then
                 do m2 = 1, 2 * Hubbard_l(nt) + 1 
                    temp = temp + v%ns( m1, m2, current_spin, na) * &
                                    proj%r(oatwfc(na)+m2, ibnd)
                 enddo
                 call daxpy (2*np, temp, swfcatom(1,oatwfc(na)+m1), 1, &
                                    hpsi(1,ibnd),              1)
              else
                 do m2 = 1, 2 * Hubbard_l(nt) + 1 
                    temp = temp + v%ns( m1, m2, current_spin, na) * &
                                    proj%k(oatwfc(na)+m2, ibnd)
                 enddo
                 call zaxpy (np, temp, swfcatom(1,oatwfc(na)+m1), 1, &
                                    hpsi(1,ibnd),              1)
              endif
           enddo
        endif
     enddo
  enddo
  call deallocate_bec_type (proj)
  !
  CALL stop_clock('vhpsi')

  return

end subroutine vhpsi


subroutine vhpsi_nc (ldap, np, mps, psip, hpsi)
  !-----------------------------------------------------------------------
  !
  ! Noncollinear version (A. Smogunov). 
  !
  USE kinds,            ONLY : DP
  USE ldaU,             ONLY : Hubbard_lmax, Hubbard_l, HUbbard_U, swfcatom, oatwfc  
  USE scf,              ONLY : v
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp
  USE noncollin_module, ONLY : npol
  USE basis,            ONLY : natomwfc
  USE wvfct,            ONLY : npwx
  USE mp_global,        ONLY : intra_bgrp_comm
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
  ALLOCATE( proj(natomwfc, mps) )

!--
! calculate <psi_at | phi_k> 
  DO ibnd = 1, mps
    DO na = 1, natomwfc
      proj(na, ibnd) = zdotc (ldap*npol, swfcatom(1, na), 1, psip(1, ibnd), 1)
    ENDDO
  ENDDO
#ifdef __MPI
  CALL mp_sum ( proj, intra_bgrp_comm )
#endif
!--

  do ibnd = 1, mps  
    do na = 1, nat  
       nt = ityp (na)  
       if (Hubbard_U(nt).ne.0.d0) then  
          nwfc = 2 * Hubbard_l(nt) + 1

          do is1 = 1, npol
           do m1 = 1, nwfc 
             temp = 0.d0
             do is2 = 1, npol
              do m2 = 1, nwfc  
                temp = temp + v%ns_nc( m1, m2, npol*(is1-1)+is2, na) * &
                              proj(oatwfc(na)+(is2-1)*nwfc+m2, ibnd)
              enddo
             enddo
             call zaxpy (ldap*npol, temp, swfcatom(1,oatwfc(na)+(is1-1)*nwfc+m1), 1, &
                         hpsi(1,ibnd),1)
           enddo
          enddo

       endif
    enddo
  enddo

  deallocate (proj)
  CALL stop_clock('vhpsi')

  return
end subroutine vhpsi_nc

